"""
Report generator for TRNspot
============================

Main class and functions for generating HTML and PDF reports.
"""

import os
from datetime import datetime
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field
from pathlib import Path
import base64

from .. import config


@dataclass
class ReportSection:
    """A section in the report."""

    title: str
    content: str
    figures: List[str] = field(default_factory=list)
    tables: List[Dict[str, Any]] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ReportData:
    """Container for all report data."""

    title: str = "TRNspot Analysis Report"
    subtitle: str = ""
    timestamp: str = field(
        default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    sections: List[ReportSection] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


class ReportGenerator:
    """
    Generate comprehensive analysis reports in HTML and PDF formats.

    Parameters
    ----------
    output_dir : str
        Directory for output files
    title : str
        Report title

    Examples
    --------
    >>> generator = ReportGenerator(output_dir="results", title="My Analysis")
    >>> generator.add_section(create_qc_section(adata))
    >>> generator.generate_html("report.html")
    >>> generator.generate_pdf("report.pdf")
    """

    def __init__(
        self,
        output_dir: str,
        title: str = "TRNspot Analysis Report",
        subtitle: str = "",
    ):
        self.output_dir = output_dir
        self.report_data = ReportData(
            title=title,
            subtitle=subtitle,
        )
        self._figures_dir = os.path.join(output_dir, "figures")

    def add_section(self, section: ReportSection) -> None:
        """Add a section to the report."""
        self.report_data.sections.append(section)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the report."""
        self.report_data.metadata[key] = value

    def _image_to_base64(self, image_path: str) -> str:
        """Convert image to base64 for embedding in HTML."""
        if not os.path.exists(image_path):
            return ""

        with open(image_path, "rb") as f:
            data = f.read()

        ext = Path(image_path).suffix.lower()
        mime_types = {
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".svg": "image/svg+xml",
            ".gif": "image/gif",
        }
        mime_type = mime_types.get(ext, "image/png")

        return f"data:{mime_type};base64,{base64.b64encode(data).decode()}"

    def _get_css(self) -> str:
        """Get CSS styles for the HTML report."""
        return """
        <style>
            :root {
                --primary-color: #2563eb;
                --secondary-color: #64748b;
                --background-color: #f8fafc;
                --card-background: #ffffff;
                --text-color: #1e293b;
                --border-color: #e2e8f0;
            }
            
            * {
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }
            
            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
                line-height: 1.6;
                color: var(--text-color);
                background-color: var(--background-color);
            }
            
            .container {
                max-width: 1200px;
                margin: 0 auto;
                padding: 2rem;
            }
            
            header {
                background: linear-gradient(135deg, var(--primary-color), #1d4ed8);
                color: white;
                padding: 3rem 2rem;
                margin-bottom: 2rem;
                border-radius: 0.5rem;
            }
            
            header h1 {
                font-size: 2.5rem;
                margin-bottom: 0.5rem;
            }
            
            header .subtitle {
                font-size: 1.2rem;
                opacity: 0.9;
            }
            
            header .timestamp {
                font-size: 0.9rem;
                opacity: 0.7;
                margin-top: 1rem;
            }
            
            .section {
                background: var(--card-background);
                border-radius: 0.5rem;
                padding: 2rem;
                margin-bottom: 2rem;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            }
            
            .section h2 {
                color: var(--primary-color);
                font-size: 1.5rem;
                margin-bottom: 1rem;
                padding-bottom: 0.5rem;
                border-bottom: 2px solid var(--border-color);
            }
            
            .section h3 {
                font-size: 1.2rem;
                margin: 1.5rem 0 0.75rem;
                color: var(--secondary-color);
            }
            
            .metrics-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 1rem;
                margin: 1rem 0;
            }
            
            .metric-card {
                background: var(--background-color);
                padding: 1rem;
                border-radius: 0.375rem;
                text-align: center;
            }
            
            .metric-card .value {
                font-size: 1.75rem;
                font-weight: 600;
                color: var(--primary-color);
            }
            
            .metric-card .label {
                font-size: 0.875rem;
                color: var(--secondary-color);
            }
            
            .figure-container {
                margin: 1.5rem 0;
                text-align: center;
            }
            
            .figure-container img {
                max-width: 100%;
                height: auto;
                border-radius: 0.375rem;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            
            .figure-container .caption {
                font-size: 0.875rem;
                color: var(--secondary-color);
                margin-top: 0.5rem;
            }
            
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 1rem 0;
            }
            
            th, td {
                padding: 0.75rem;
                text-align: left;
                border-bottom: 1px solid var(--border-color);
            }
            
            th {
                background: var(--background-color);
                font-weight: 600;
            }
            
            tr:hover {
                background-color: var(--background-color);
            }
            
            .toc {
                background: var(--card-background);
                padding: 1.5rem;
                border-radius: 0.5rem;
                margin-bottom: 2rem;
            }
            
            .toc h3 {
                margin-bottom: 1rem;
            }
            
            .toc ul {
                list-style: none;
            }
            
            .toc li {
                margin: 0.5rem 0;
            }
            
            .toc a {
                color: var(--primary-color);
                text-decoration: none;
            }
            
            .toc a:hover {
                text-decoration: underline;
            }
            
            footer {
                text-align: center;
                padding: 2rem;
                color: var(--secondary-color);
                font-size: 0.875rem;
            }
            
            @media print {
                body { background: white; }
                .container { max-width: none; padding: 0; }
                .section { box-shadow: none; page-break-inside: avoid; }
                header { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
            }
        </style>
        """

    def _render_metrics(self, metrics: Dict[str, Any]) -> str:
        """Render metrics as HTML cards."""
        if not metrics:
            return ""

        cards = []
        for label, value in metrics.items():
            if isinstance(value, float):
                value_str = f"{value:.2f}"
            else:
                value_str = str(value)

            cards.append(
                f"""
            <div class="metric-card">
                <div class="value">{value_str}</div>
                <div class="label">{label}</div>
            </div>
            """
            )

        return f'<div class="metrics-grid">{"".join(cards)}</div>'

    def _render_figures(self, figures: List[str], embed: bool = True) -> str:
        """Render figures as HTML."""
        if not figures:
            return ""

        html_parts = []
        for fig_path in figures:
            if not os.path.exists(fig_path):
                continue

            if embed:
                src = self._image_to_base64(fig_path)
            else:
                src = fig_path

            caption = Path(fig_path).stem.replace("_", " ").title()
            html_parts.append(
                f"""
            <div class="figure-container">
                <img src="{src}" alt="{caption}">
                <div class="caption">{caption}</div>
            </div>
            """
            )

        return "".join(html_parts)

    def _render_tables(self, tables: List[Dict[str, Any]]) -> str:
        """Render tables as HTML."""
        if not tables:
            return ""

        html_parts = []
        for table_data in tables:
            title = table_data.get("title", "")
            headers = table_data.get("headers", [])
            rows = table_data.get("rows", [])

            if title:
                html_parts.append(f"<h3>{title}</h3>")

            header_html = "".join(f"<th>{h}</th>" for h in headers)
            rows_html = ""
            for row in rows:
                cells = "".join(f"<td>{cell}</td>" for cell in row)
                rows_html += f"<tr>{cells}</tr>"

            html_parts.append(
                f"""
            <table>
                <thead><tr>{header_html}</tr></thead>
                <tbody>{rows_html}</tbody>
            </table>
            """
            )

        return "".join(html_parts)

    def _render_section(self, section: ReportSection, idx: int) -> str:
        """Render a single section as HTML."""
        section_id = f"section-{idx}"

        return f"""
        <div class="section" id="{section_id}">
            <h2>{section.title}</h2>
            {self._render_metrics(section.metrics)}
            <div class="content">{section.content}</div>
            {self._render_figures(section.figures)}
            {self._render_tables(section.tables)}
        </div>
        """

    def _render_toc(self) -> str:
        """Render table of contents."""
        items = []
        for idx, section in enumerate(self.report_data.sections):
            items.append(f'<li><a href="#section-{idx}">{section.title}</a></li>')

        return f"""
        <div class="toc">
            <h3>Table of Contents</h3>
            <ul>{"".join(items)}</ul>
        </div>
        """

    def generate_html(
        self,
        output_path: Optional[str] = None,
        embed_images: bool = True,
    ) -> str:
        """
        Generate HTML report.

        Parameters
        ----------
        output_path : str, optional
            Path to save the HTML file
        embed_images : bool
            Whether to embed images as base64

        Returns
        -------
        str
            HTML content
        """
        sections_html = "\n".join(
            self._render_section(section, idx)
            for idx, section in enumerate(self.report_data.sections)
        )

        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{self.report_data.title}</title>
            {self._get_css()}
        </head>
        <body>
            <div class="container">
                <header>
                    <h1>{self.report_data.title}</h1>
                    <div class="subtitle">{self.report_data.subtitle}</div>
                    <div class="timestamp">Generated: {self.report_data.timestamp}</div>
                </header>
                
                {self._render_toc()}
                
                {sections_html}
                
                <footer>
                    Generated by TRNspot v{config.__version__ if hasattr(config, '__version__') else '0.1.0'}
                </footer>
            </div>
        </body>
        </html>
        """

        if output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html)
            print(f"✓ HTML report saved to: {output_path}")

        return html

    def generate_pdf(
        self,
        output_path: Optional[str] = None,
    ) -> Optional[str]:
        """
        Generate PDF report from HTML.

        Requires weasyprint or pdfkit to be installed.

        Parameters
        ----------
        output_path : str, optional
            Path to save the PDF file

        Returns
        -------
        str or None
            Path to PDF file, or None if generation failed
        """
        if output_path is None:
            output_path = os.path.join(self.output_dir, "report.pdf")

        html_content = self.generate_html(embed_images=True)

        # Try weasyprint first
        try:
            from weasyprint import HTML

            HTML(string=html_content).write_pdf(output_path)
            print(f"✓ PDF report saved to: {output_path}")
            return output_path
        except ImportError:
            pass

        # Try pdfkit as fallback
        try:
            import pdfkit

            pdfkit.from_string(html_content, output_path)
            print(f"✓ PDF report saved to: {output_path}")
            return output_path
        except ImportError:
            pass

        print("⚠ PDF generation requires 'weasyprint' or 'pdfkit'. Install with:")
        print("  pip install weasyprint")
        print("  or")
        print("  pip install pdfkit")
        return None


def generate_report(
    output_dir: str,
    title: str = "TRNspot Analysis Report",
    adata=None,
    celloracle_result=None,
    hotspot_result=None,
    formats: List[str] = ["html", "pdf"],
) -> Dict[str, str]:
    """
    Generate analysis report in specified formats.

    Parameters
    ----------
    output_dir : str
        Output directory
    title : str
        Report title
    adata : AnnData, optional
        Processed AnnData object
    celloracle_result : tuple, optional
        CellOracle results (oracle, links)
    hotspot_result : Hotspot, optional
        Hotspot results
    formats : list
        Output formats ('html', 'pdf')

    Returns
    -------
    dict
        Paths to generated report files
    """
    from .sections import (
        create_qc_section,
        create_preprocessing_section,
        create_celloracle_section,
        create_hotspot_section,
        create_summary_section,
    )

    generator = ReportGenerator(output_dir=output_dir, title=title)

    # Add sections based on available data
    if adata is not None:
        generator.add_section(create_qc_section(adata, output_dir))
        generator.add_section(create_preprocessing_section(adata, output_dir))

    if celloracle_result is not None:
        generator.add_section(create_celloracle_section(celloracle_result, output_dir))

    if hotspot_result is not None:
        generator.add_section(create_hotspot_section(hotspot_result, output_dir))

    generator.add_section(
        create_summary_section(
            adata=adata,
            celloracle_result=celloracle_result,
            hotspot_result=hotspot_result,
            output_dir=output_dir,
        )
    )

    # Generate reports
    outputs = {}

    if "html" in formats:
        html_path = os.path.join(output_dir, "report.html")
        generator.generate_html(html_path)
        outputs["html"] = html_path

    if "pdf" in formats:
        pdf_path = os.path.join(output_dir, "report.pdf")
        result = generator.generate_pdf(pdf_path)
        if result:
            outputs["pdf"] = pdf_path

    return outputs


def generate_html_report(
    output_dir: str,
    **kwargs,
) -> str:
    """Convenience function to generate HTML report only."""
    outputs = generate_report(output_dir, formats=["html"], **kwargs)
    return outputs.get("html", "")


def generate_pdf_report(
    output_dir: str,
    **kwargs,
) -> Optional[str]:
    """Convenience function to generate PDF report only."""
    outputs = generate_report(output_dir, formats=["pdf"], **kwargs)
    return outputs.get("pdf")
