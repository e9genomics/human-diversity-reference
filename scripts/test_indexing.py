from remap_divref import Haplotype


def create_haplotype(
    sequence_id="test_hap",
    sequence="ACGT",
    sequence_length=4,
    n_variants=3,
    variants="1:500:A:T,1:505:C:G,1:510:T:A",
    gnomAD_AF_afr="0.1,0.2,0.3",
    gnomAD_AF_amr="0.15,0.25,0.35",
    gnomAD_AF_eas="0.12,0.22,0.32",
    gnomAD_AF_nfe="0.11,0.21,0.31",
    gnomAD_AF_sas="0.13,0.23,0.33",
    **kwargs,
):
    defaults = {
        "fraction_phased": 1.0,
        "popmax_empirical_AF": 0.25,
        "popmax_empirical_AN": 1000,
        "estimated_gnomad_AF": 0.15,
        "source": "test_source",
    }
    # Update defaults with any provided kwargs
    defaults.update(kwargs)

    return Haplotype(
        sequence_id=sequence_id,
        sequence=sequence,
        sequence_length=sequence_length,
        n_variants=n_variants,
        variants=variants,
        gnomAD_AF_afr=gnomAD_AF_afr,
        gnomAD_AF_amr=gnomAD_AF_amr,
        gnomAD_AF_eas=gnomAD_AF_eas,
        gnomAD_AF_nfe=gnomAD_AF_nfe,
        gnomAD_AF_sas=gnomAD_AF_sas,
        **defaults,
    )


def test_variants_involved():
    # Test case 1: Basic single-base variants
    haplotype1 = create_haplotype(
        sequence_id="hap1", variants="1:500:A:T,1:505:C:G,1:510:T:A"
    )
    # With context_size=10, pos 0 in sequence = pos 490 in GRCh38
    # Looking at sequence positions 12-17 (GRCh38 502-507)
    ref_mapping = haplotype1.reference_mapping(12, 17, 10)
    assert ref_mapping.first_variant_index == 1
    assert ref_mapping.last_variant_index == 1  # Should only include the second variant
    assert ref_mapping.gnomad_popmax_pop == "amr"  # AMR has highest minimum frequency
    assert ref_mapping.gnomad_popmax_frequency == 0.15  # Minimum AMR frequency

    # Check population frequencies
    expected_pop_freqs = {
        "afr": [0.2],
        "amr": [0.25],
        "eas": [0.22],
        "nfe": [0.21],
        "sas": [0.23],
    }
    assert ref_mapping.population_frequencies == expected_pop_freqs

    # Test case 2: Deletion affects coordinate mapping
    haplotype2 = create_haplotype(
        sequence_id="hap2", variants="1:500:A:T,1:505:CC:G,1:510:T:A"
    )
    ref_mapping = haplotype2.reference_mapping(14, 20, 10)
    assert ref_mapping.first_variant_index == 1
    assert (
        ref_mapping.last_variant_index == 2
    )  # Should include second and third variants
    assert ref_mapping.gnomad_popmax_pop == "amr"
    assert ref_mapping.gnomad_popmax_frequency == 0.15

    # Test case 3: Insertion affects coordinate mapping
    haplotype3 = create_haplotype(
        sequence_id="hap3", variants="1:500:T:TTT,1:505:C:G,1:510:T:A"
    )
    ref_mapping = haplotype3.reference_mapping(15, 18, 10)
    assert ref_mapping.first_variant_index == 1
    assert ref_mapping.last_variant_index == 1  # Should only include second variant
    assert ref_mapping.gnomad_popmax_pop == "amr"
    assert ref_mapping.gnomad_popmax_frequency == 0.15

    # Test case 4: No variants in range
    haplotype4 = create_haplotype(sequence_id="hap4", variants="1:500:A:T,1:510:C:G")
    ref_mapping = haplotype4.reference_mapping(0, 5, 10)
    assert ref_mapping.first_variant_index is None
    assert ref_mapping.last_variant_index is None
    assert ref_mapping.gnomad_popmax_pop == "amr"
    assert ref_mapping.gnomad_popmax_frequency == 0.15

    # Test case 5: Multiple indels with complex shifts and different population frequencies
    haplotype5 = create_haplotype(
        sequence_id="hap5",
        variants="1:500:AT:A,1:505:C:CTT,1:510:GGG:T",
        gnomAD_AF_afr="0.05,0.15,0.25",  # min = 0.05
        gnomAD_AF_amr="0.06,0.16,0.26",  # min = 0.06
        gnomAD_AF_eas="0.07,0.17,0.27",  # min = 0.07 (highest min)
        gnomAD_AF_nfe="0.04,0.14,0.24",  # min = 0.04
        gnomAD_AF_sas="0.03,0.13,0.23",  # min = 0.03
    )
    ref_mapping = haplotype5.reference_mapping(9, 22, 10)
    assert ref_mapping.first_variant_index == 0
    assert ref_mapping.last_variant_index == 2  # Should include all variants
    assert ref_mapping.gnomad_popmax_pop == "eas"  # EAS has highest minimum frequency
    assert ref_mapping.gnomad_popmax_frequency == 0.07  # Minimum EAS frequency

    # Test case 6: specific insertion with null frequencies
    haplotype6 = create_haplotype(
        sequence_id="hap7",
        sequence="ACTACTATCTATATCATCTACTACTACTATCATCATCATCAT",
        sequence_length=4,
        n_variants=1,
        variants="chr12:90349349:T:TATGCAAGTGTCATCAGATGAATTGATGACATTTTTGTCAAGTTTAAGCACTGAAAGAACAAACCTCTAAATC",
        gnomAD_AF_afr="null",
        gnomAD_AF_amr="null",
        gnomAD_AF_eas="null",
        gnomAD_AF_nfe="null",
        gnomAD_AF_sas="null",
    )
    ref_mapping = haplotype6.reference_mapping(22, 42, 25)
    assert ref_mapping.first_variant_index == 0
    assert ref_mapping.last_variant_index == 0
    assert ref_mapping.start == 90349346
    assert ref_mapping.end == 90349350
    assert ref_mapping.gnomad_popmax_frequency == 0

    # Population frequencies should now be a dictionary with lists
    expected_null_freqs = {"afr": [0], "amr": [0], "eas": [0], "nfe": [0], "sas": [0]}
    assert ref_mapping.population_frequencies == expected_null_freqs
