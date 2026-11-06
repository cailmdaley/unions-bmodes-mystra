"""Regression checks for report navigation without building the publication."""
import unittest

from check_report_links import broken_links
from check_report_build import check


class ReportLinksTest(unittest.TestCase):
    def test_myst_and_equation_labels(self):
        text = r"""
(sec:results)=
(data-availability .unnumbered)=
$$ (eq:transform)
\label{eq:estimator}
[results](#sec:results) [data](#data-availability)
[transform](#eq:transform) [estimator](#eq:estimator)
[inventory](#astra-inventory)
"""
        self.assertEqual(broken_links(text), [])

    def test_missing_sections_and_raw_record_fragments_fail(self):
        self.assertEqual(
            broken_links("[output](#outputs.pte_table_appendix) [missing](#sec:missing)"),
            ["outputs.pte_table_appendix", "sec:missing"],
        )


class ReportBuildTest(unittest.TestCase):
    def test_bound_value_survives_build(self):
        source = '{astra:value col=display where="macro=x"}`outputs.paper_claims`'
        ast = {"type": "span", "data": {"astra": {"kind": "value"}}}
        self.assertEqual(check(source, ast), [])

    def test_caption_role_left_unresolved_fails(self):
        source = '{astra:value col=display where="macro=x"}`outputs.paper_claims`'
        ast = {"type": "mystRole", "name": "astra:value"}
        self.assertEqual(check(source, ast), [
            "bound values: 0 resolved, 1 authored", "unresolved ASTRA roles: 1"
        ])


if __name__ == "__main__":
    unittest.main()
