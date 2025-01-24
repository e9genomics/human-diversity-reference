from unittest.mock import patch

import hail as hl

from create_fasta_and_index import get_haplo_sequence


# Helper function to create test variants
def create_variant(contig: str, position: int, ref: str, alt: str) -> hl.Struct:
    """
    Creates a test variant with the specified paramet ers.
    Returns a struct matching the format expected by the haplotype functions.
    """
    return hl.Struct(
        locus=hl.Struct(contig=contig, position=position),
        alleles=[ref, alt]
    )


def create_reference_mock(reference_sequence):
    """
    Creates a mock function that behaves like hl.get_sequence but uses a predefined reference.

    The function takes a reference string and returns a mock that can handle genomic coordinates
    by returning the appropriate subsequence, just like a real reference genome would.
    """

    def mock_get_sequence(contig, position, before=0, after=0, reference_genome=None):
        """
        Mock implementation of hl.get_sequence that returns subsequences from our reference.

        Args:
            contig: Chromosome (ignored in this mock)
            position: 1-based position in the reference
            before: Number of bases to include before the position
            after: Number of bases to include after the position
        """
        # Convert to 0-based indexing
        return hl.str(reference_sequence)[position - before:position + after + 1]

    return mock_get_sequence


def test_get_haplo_sequence_edge_cases():
    reference = "01234567891"

    two_snps = [
        create_variant("chr1", 4, "A", "T"),
        create_variant("chr1", 6, "G", "C")
    ]

    two_insertions = [
        create_variant("chr1", 4, "A", "AT"),
        create_variant("chr1", 6, "G", "GC"),
    ]

    two_deletions = [
        create_variant("chr1", 4, "AT", "A"),
        create_variant("chr1", 7, "GC", "G")
    ]

    mock_get_sequence = create_reference_mock(reference)

    with patch('hail.get_sequence', side_effect=mock_get_sequence):
        result_snps = hl.eval(get_haplo_sequence(context_size=2, variants=two_snps))
        assert result_snps == '23T5C78'

        result_insertions = hl.eval(get_haplo_sequence(context_size=2, variants=two_insertions))
        assert result_insertions == '23AT5GC78'

        result_deletions = hl.eval(get_haplo_sequence(context_size=2, variants=two_deletions))
        assert result_deletions == '23A6G91'
