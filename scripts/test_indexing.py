import pytest
from remap_divref import Haplotype

def test_variants_involved():
    # Test case 1: Basic single-base variants
    haplotype1 = Haplotype(
        sequence_id="hap1",
        sequence="ACGT",
        sequence_length=4,
        n_variants=3,
        variants="1:500:A:T,1:505:C:G,1:510:T:A",
        gnomad_AF="0.1,0.2,0.3"
    )
    # With context_size=10, pos 0 in sequence = pos 490 in GRCh38
    # Looking at sequence positions 12-17 (GRCh38 502-507)
    ref_mapping = haplotype1.reference_mapping(12, 17, 10)
    assert ref_mapping.first_variant_index == 1
    assert ref_mapping.last_variant_index == 1  # Should only include the second variant

    # Test case 2: Deletion affects coordinate mapping
    haplotype2 = Haplotype(
        sequence_id="hap2",
        sequence="ACGT",
        sequence_length=4,
        n_variants=3,
        variants="1:500:A:T,1:505:CC:G,1:510:T:A",
        gnomad_AF="0.1,0.2,0.3"
    )
    # With context_size=10:
    # First variant at seq pos 10 (no shift)
    # Second variant at seq pos 15 (2bp becomes 1bp, -1 shift)
    # Third variant at seq pos 19 (shifted left by 1 from previous deletion)
    ref_mapping = haplotype2.reference_mapping(14, 20, 10)
    assert ref_mapping.first_variant_index == 1
    assert ref_mapping.last_variant_index == 2  # Should include second and third variants

    # Test case 3: Insertion affects coordinate mapping
    haplotype3 = Haplotype(
        sequence_id="hap3",
        sequence="ACGT",
        sequence_length=4,
        n_variants=3,
        variants="1:500:T:TTT,1:505:C:G,1:510:T:A",
        gnomad_AF="0.1,0.2,0.3"
    )
    # With context_size=10:
    # First variant at seq pos 10 (1bp becomes 3bp, +2 shift)
    # Second variant at seq pos 17 (shifted right by 2)
    # Third variant at seq pos 22 (shifted right by 2)
    ref_mapping = haplotype3.reference_mapping(15, 18, 10)
    assert ref_mapping.first_variant_index == 1
    assert ref_mapping.last_variant_index == 1  # Should only include second variant

    # Test case 4: No variants in range
    haplotype4 = Haplotype(
        sequence_id="hap4",
        sequence="ACGT",
        sequence_length=4,
        n_variants=2,
        variants="1:500:A:T,1:510:C:G",
        gnomad_AF="0.1,0.2"
    )
    # Looking at sequence range that doesn't overlap any variants
    ref_mapping = haplotype4.reference_mapping(0, 5, 10)
    assert ref_mapping.first_variant_index == None
    assert ref_mapping.last_variant_index == None

    # Test case 5: Multiple indels with complex shifts
    haplotype5 = Haplotype(
        sequence_id="hap5",
        sequence="ACGT",
        sequence_length=4,
        n_variants=3,
        variants="1:500:AT:A,1:505:C:CTT,1:510:GGG:T",
        gnomad_AF="0.1,0.2,0.3"
    )
    # With context_size=10:
    # First variant at seq pos 10 (2bp becomes 1bp, -1 shift)
    # Second variant at seq pos 14 (1bp becomes 3bp, +2 shift)
    # Third variant at seq pos 21 (shifted right by 1 total, then 3bp becomes 1bp)
    ref_mapping = haplotype5.reference_mapping(9, 22, 10)
    assert ref_mapping.first_variant_index == 0
    assert ref_mapping.last_variant_index == 2  # Should include all variants

    # Test case 6: specific insertion
    haplotype6 = Haplotype(
        sequence_id="hap7",  # unused
        sequence="ACTACTATCTATATCATCTACTACTACTATCATCATCATCAT",  # unused
        sequence_length=4,  # unused
        n_variants=1,
        variants="chr12:90349349:T:TATGCAAGTGTCATCAGATGAATTGATGACATTTTTGTCAAGTTTAAGCACTGAAAGAACAAACCTCTAAATC",
        gnomad_AF="0.1,0.2,0.3"  # unused
    )
    # Looking at sequence range that doesn't overlap any variants
    ref_mapping = haplotype6.reference_mapping(22, 42, 25)
    assert ref_mapping.first_variant_index == 0
    assert ref_mapping.last_variant_index == 0
    assert ref_mapping.start == 90349346
    assert ref_mapping.end == 90349350
