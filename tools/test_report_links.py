"""Regression checks for report navigation without building the publication."""
import unittest

from check_report_links import broken_links


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
            broken_links("[output](#outputs.selected_config) [missing](#sec:missing)"),
            ["outputs.selected_config", "sec:missing"],
        )


if __name__ == "__main__":
    unittest.main()
