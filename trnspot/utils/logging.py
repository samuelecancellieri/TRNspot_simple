"""
Logging utility functions for TRNspot
=====================================

Functions for pipeline logging and error tracking.
"""

import os
import logging
from typing import Optional, Dict, Any, Tuple

# Global logger instances
_pipeline_logger: Optional[logging.Logger] = None
_error_logger: Optional[logging.Logger] = None


def setup_logging(
    output_dir: str,
    pipeline_log: str = "pipeline.log",
    error_log: str = "error.log",
) -> Tuple[logging.Logger, logging.Logger]:
    """
    Setup logging system for pipeline execution and errors.

    Creates two log files:
    - pipeline.log: Records all pipeline steps with timestamps
    - error.log: Records all errors with full tracebacks

    Parameters
    ----------
    output_dir : str
        Directory where log files will be created
    pipeline_log : str
        Name of the pipeline log file
    error_log : str
        Name of the error log file

    Returns
    -------
    tuple
        (pipeline_logger, error_logger)
    """
    global _pipeline_logger, _error_logger

    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Pipeline logger - tracks all steps
    _pipeline_logger = logging.getLogger("trnspot.pipeline")
    _pipeline_logger.setLevel(logging.INFO)
    _pipeline_logger.handlers.clear()

    pipeline_handler = logging.FileHandler(
        os.path.join(log_dir, pipeline_log), mode="a"
    )
    pipeline_handler.setLevel(logging.INFO)
    pipeline_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    pipeline_handler.setFormatter(pipeline_formatter)
    _pipeline_logger.addHandler(pipeline_handler)

    # Error logger - tracks all errors with tracebacks
    _error_logger = logging.getLogger("trnspot.error")
    _error_logger.setLevel(logging.ERROR)
    _error_logger.handlers.clear()

    error_handler = logging.FileHandler(os.path.join(log_dir, error_log), mode="a")
    error_handler.setLevel(logging.ERROR)
    error_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s\n%(exc_info)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    error_handler.setFormatter(error_formatter)
    _error_logger.addHandler(error_handler)

    # Log session start
    _pipeline_logger.info("=" * 70)
    _pipeline_logger.info("NEW PIPELINE SESSION STARTED")
    _pipeline_logger.info("=" * 70)

    return _pipeline_logger, _error_logger


def get_pipeline_logger() -> Optional[logging.Logger]:
    """Get the global pipeline logger instance."""
    return _pipeline_logger


def get_error_logger() -> Optional[logging.Logger]:
    """Get the global error logger instance."""
    return _error_logger


def log_step(
    step_name: str,
    status: str = "STARTED",
    details: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Log a pipeline step with timestamp.

    Parameters
    ----------
    step_name : str
        Name of the pipeline step
    status : str
        Status of the step (STARTED, COMPLETED, FAILED, SKIPPED, etc.)
    details : dict, optional
        Additional details to log (e.g., n_obs, n_vars)

    Examples
    --------
    >>> log_step("Preprocessing", "STARTED", {"n_cells": 1000})
    >>> log_step("Preprocessing", "COMPLETED")
    """
    if _pipeline_logger is None:
        return

    message = f"[{step_name}] {status}"
    if details:
        detail_str = ", ".join([f"{k}={v}" for k, v in details.items()])
        message += f" - {detail_str}"

    _pipeline_logger.info(message)


def log_error(
    error_context: str,
    exception: Exception,
) -> None:
    """
    Log an error with full traceback.

    Parameters
    ----------
    error_context : str
        Description of where/when the error occurred
    exception : Exception
        The exception that was raised

    Examples
    --------
    >>> try:
    ...     risky_operation()
    ... except Exception as e:
    ...     log_error("RiskyOperation", e)
    ...     raise
    """
    if _error_logger is None:
        return

    _error_logger.error(
        f"ERROR in {error_context}: {str(exception)}",
        exc_info=True,
    )

    # Also log to pipeline logger
    if _pipeline_logger:
        _pipeline_logger.error(f"ERROR in {error_context}: {str(exception)}")
