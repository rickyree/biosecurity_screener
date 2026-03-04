#!/usr/bin/env python3
"""
bioRxiv Preprint Client

Python interface for querying the bioRxiv preprint repository.
Supports searching by keywords, authors, date ranges, and subject categories.
"""

import argparse
import json
import sys
import time
from datetime import datetime, timedelta
from typing import Optional

import requests


# All valid bioRxiv subject categories
SUBJECT_AREAS = (
    "animal-behavior-and-cognition",
    "biochemistry",
    "bioengineering",
    "bioinformatics",
    "biophysics",
    "cancer-biology",
    "cell-biology",
    "clinical-trials",
    "developmental-biology",
    "ecology",
    "epidemiology",
    "evolutionary-biology",
    "genetics",
    "genomics",
    "immunology",
    "microbiology",
    "molecular-biology",
    "neuroscience",
    "paleontology",
    "pathology",
    "pharmacology-and-toxicology",
    "physiology",
    "plant-biology",
    "scientific-communication-and-education",
    "synthetic-biology",
    "systems-biology",
    "zoology",
)


class PreprintClient:
    """Client for accessing the bioRxiv preprint API."""

    API_ROOT = "https://api.biorxiv.org"
    REQUEST_DELAY = 0.5  # Seconds between requests

    def __init__(self, debug: bool = False):
        self._debug = debug
        self._http = requests.Session()
        self._http.headers["User-Agent"] = "PreprintClient/1.0"

    def _print_debug(self, msg: str) -> None:
        if self._debug:
            print(f"[DEBUG] {msg}", file=sys.stderr)

    def _request(self, path: str) -> dict:
        """Execute an API request with rate limiting."""
        url = f"{self.API_ROOT}/{path}"
        self._print_debug(f"GET {url}")

        try:
            resp = self._http.get(url, timeout=30)
            resp.raise_for_status()
            time.sleep(self.REQUEST_DELAY)
            return resp.json()
        except requests.RequestException as err:
            self._print_debug(f"Request failed: {err}")
            return {"messages": [{"status": "error"}], "collection": []}

    def find_by_date(
        self,
        since: str,
        until: str,
        subject: Optional[str] = None,
    ) -> list:
        """
        Retrieve preprints posted within a date range.

        Args:
            since: Start date (YYYY-MM-DD)
            until: End date (YYYY-MM-DD)
            subject: Optional subject category

        Returns:
            List of preprint dictionaries
        """
        self._print_debug(f"Date search: {since} to {until}")

        if subject and subject in SUBJECT_AREAS:
            path = f"details/biorxiv/{since}/{until}/{subject}"
        else:
            path = f"details/biorxiv/{since}/{until}"

        result = self._request(path)
        papers = result.get("collection", [])
        self._print_debug(f"Retrieved {len(papers)} preprints")
        return papers

    def find_by_author(
        self,
        name: str,
        since: Optional[str] = None,
        until: Optional[str] = None,
    ) -> list:
        """
        Find preprints by author name.

        Args:
            name: Author name (partial match, case-insensitive)
            since: Start date (defaults to 3 years ago)
            until: End date (defaults to today)

        Returns:
            List of matching preprints
        """
        # Default to 3-year window
        if until is None:
            until = datetime.now().strftime("%Y-%m-%d")
        if since is None:
            since = (datetime.now() - timedelta(days=3 * 365)).strftime("%Y-%m-%d")

        self._print_debug(f"Author search: '{name}'")

        all_papers = self.find_by_date(since, until)
        name_lower = name.lower()

        matches = [
            p for p in all_papers
            if name_lower in p.get("authors", "").lower()
        ]

        self._print_debug(f"Found {len(matches)} papers by '{name}'")
        return matches

    def find_by_terms(
        self,
        terms: list,
        since: Optional[str] = None,
        until: Optional[str] = None,
        subject: Optional[str] = None,
        fields: Optional[list] = None,
    ) -> list:
        """
        Search for preprints containing keywords.

        Args:
            terms: Keywords to search for
            since: Start date (defaults to 1 year ago)
            until: End date (defaults to today)
            subject: Optional subject filter
            fields: Which fields to search (default: title, abstract)

        Returns:
            List of matching preprints
        """
        # Default to 1-year window
        if until is None:
            until = datetime.now().strftime("%Y-%m-%d")
        if since is None:
            since = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d")

        if fields is None:
            fields = ["title", "abstract"]

        self._print_debug(f"Term search: {terms}")

        all_papers = self.find_by_date(since, until, subject)
        terms_lower = [t.lower() for t in terms]

        matches = []
        for paper in all_papers:
            searchable = " ".join(
                str(paper.get(f, "")).lower() for f in fields
            )
            if any(term in searchable for term in terms_lower):
                matches.append(paper)

        self._print_debug(f"Found {len(matches)} papers matching terms")
        return matches

    def get_by_doi(self, doi: str) -> dict:
        """
        Fetch metadata for a specific preprint.

        Args:
            doi: The preprint DOI (with or without URL prefix)

        Returns:
            Preprint metadata dictionary, or empty dict if not found
        """
        # Strip URL prefix if present
        if "doi.org/" in doi:
            doi = doi.split("doi.org/")[-1]

        self._print_debug(f"DOI lookup: {doi}")

        result = self._request(f"details/biorxiv/{doi}")
        collection = result.get("collection", [])

        return collection[0] if collection else {}

    def fetch_pdf(self, doi: str, destination: str) -> bool:
        """
        Download a preprint PDF.

        Args:
            doi: The preprint DOI
            destination: Local file path for the PDF

        Returns:
            True if successful, False otherwise
        """
        if "doi.org/" in doi:
            doi = doi.split("doi.org/")[-1]

        pdf_url = f"https://www.biorxiv.org/content/{doi}v1.full.pdf"
        self._print_debug(f"Downloading: {pdf_url}")

        try:
            resp = self._http.get(pdf_url, timeout=60)
            resp.raise_for_status()

            with open(destination, "wb") as f:
                f.write(resp.content)

            self._print_debug(f"Saved to: {destination}")
            return True

        except Exception as err:
            self._print_debug(f"Download failed: {err}")
            return False

    def normalize(self, paper: dict, include_abstract: bool = True) -> dict:
        """
        Convert API response to standardized format.

        Args:
            paper: Raw paper dict from API
            include_abstract: Whether to include abstract text

        Returns:
            Normalized dictionary with consistent field names
        """
        doi = paper.get("doi", "")
        version = paper.get("version", "1")

        normalized = {
            "doi": doi,
            "title": paper.get("title", ""),
            "authors": paper.get("authors", ""),
            "corresponding_author": paper.get("author_corresponding", ""),
            "institution": paper.get("author_corresponding_institution", ""),
            "posted": paper.get("date", ""),
            "revision": version,
            "category": paper.get("category", ""),
            "license": paper.get("license", ""),
            "paper_type": paper.get("type", ""),
            "journal_ref": paper.get("published", ""),
        }

        if include_abstract:
            normalized["abstract"] = paper.get("abstract", "")

        if doi:
            normalized["pdf_link"] = f"https://www.biorxiv.org/content/{doi}v{version}.full.pdf"
            normalized["web_link"] = f"https://www.biorxiv.org/content/{doi}v{version}"

        return normalized


def _parse_args():
    """Set up command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Query bioRxiv preprints",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug output",
    )

    # Query type
    query_grp = parser.add_argument_group("Query")
    query_grp.add_argument(
        "-t", "--terms",
        nargs="+",
        metavar="TERM",
        help="Search terms (keywords)",
    )
    query_grp.add_argument(
        "-a", "--author",
        metavar="NAME",
        help="Author name to search",
    )
    query_grp.add_argument(
        "--doi",
        metavar="DOI",
        help="Specific DOI to retrieve",
    )

    # Date range
    date_grp = parser.add_argument_group("Date Range")
    date_grp.add_argument(
        "--since",
        metavar="DATE",
        help="Start date (YYYY-MM-DD)",
    )
    date_grp.add_argument(
        "--until",
        metavar="DATE",
        help="End date (YYYY-MM-DD)",
    )
    date_grp.add_argument(
        "--recent",
        type=int,
        metavar="DAYS",
        help="Search last N days",
    )

    # Filters
    filter_grp = parser.add_argument_group("Filters")
    filter_grp.add_argument(
        "-s", "--subject",
        choices=SUBJECT_AREAS,
        metavar="CATEGORY",
        help="Subject category filter",
    )
    filter_grp.add_argument(
        "--fields",
        nargs="+",
        choices=["title", "abstract", "authors"],
        default=["title", "abstract"],
        help="Fields to search for keywords",
    )

    # Output
    out_grp = parser.add_argument_group("Output")
    out_grp.add_argument(
        "-o", "--out",
        metavar="FILE",
        help="Output file (default: stdout)",
    )
    out_grp.add_argument(
        "--max",
        type=int,
        metavar="N",
        help="Maximum results to return",
    )
    out_grp.add_argument(
        "--fetch-pdf",
        metavar="PATH",
        help="Download PDF (requires --doi)",
    )

    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = _parse_args()

    client = PreprintClient(debug=args.verbose)

    # Resolve date range
    end_date = args.until or datetime.now().strftime("%Y-%m-%d")
    if args.recent:
        start_date = (datetime.now() - timedelta(days=args.recent)).strftime("%Y-%m-%d")
    else:
        start_date = args.since

    # Handle PDF download
    if args.fetch_pdf:
        if not args.doi:
            print("Error: --doi required for PDF download", file=sys.stderr)
            return 1
        ok = client.fetch_pdf(args.doi, args.fetch_pdf)
        return 0 if ok else 1

    # Execute query
    papers = []

    if args.doi:
        paper = client.get_by_doi(args.doi)
        if paper:
            papers = [paper]

    elif args.author:
        papers = client.find_by_author(args.author, start_date, end_date)

    elif args.terms:
        if not start_date:
            print("Error: --since or --recent required for term search", file=sys.stderr)
            return 1
        papers = client.find_by_terms(
            args.terms, start_date, end_date, args.subject, args.fields
        )

    elif start_date:
        papers = client.find_by_date(start_date, end_date, args.subject)

    else:
        print("Error: Specify --terms, --author, --doi, or date range", file=sys.stderr)
        return 1

    # Apply limit
    if args.max and args.max > 0:
        papers = papers[:args.max]

    # Normalize output
    normalized = [client.normalize(p) for p in papers]

    output = {
        "query": {
            "terms": args.terms,
            "author": args.author,
            "doi": args.doi,
            "since": start_date,
            "until": end_date,
            "subject": args.subject,
        },
        "count": len(normalized),
        "papers": normalized,
    }

    output_str = json.dumps(output, indent=2)

    if args.out:
        with open(args.out, "w") as f:
            f.write(output_str)
        print(f"Wrote {len(normalized)} results to {args.out}", file=sys.stderr)
    else:
        print(output_str)

    return 0


if __name__ == "__main__":
    sys.exit(main())
