"""
TRNspot reporting module
========================

This module generates comprehensive PDF and HTML reports from pipeline results.
"""

from .generator import (
    ReportGenerator,
    generate_report,
    generate_html_report,
    generate_pdf_report,
)
from .sections import (
    create_qc_section,
    create_preprocessing_section,
    create_celloracle_section,
    create_hotspot_section,
    create_summary_section,
)

__all__ = [
    # Main generator
    "ReportGenerator",
    "generate_report",
    "generate_html_report",
    "generate_pdf_report",
    # Section builders
    "create_qc_section",
    "create_preprocessing_section",
    "create_celloracle_section",
    "create_hotspot_section",
    "create_summary_section",
]
