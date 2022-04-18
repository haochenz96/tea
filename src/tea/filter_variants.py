def filter_variants(self, min_dp=10, min_gq=30, min_vaf=20, max_vaf=100, min_prct_cells=50, min_mut_prct_cells=1, min_std=0):
    """
    Find informative variants.

    This method also adds the `NGT_FILTERED` layer to the assay
    which is a copy of the NGT layer but with the NGT for the
    cell-variants not passing the filters set to 3 i.e. missing.

    Parameters
    ----------
    min_dp : int
        The minimum depth (DP) for the call to be considered.
        Variants with less than this DP in a given
        barcode are treated as no calls.
    min_gq : int
        The minimum genotype quality (GQ) for the call to be
        considered. Variants with less than this GQ
        in a given barcode are treated as no calls.
    min_vaf : float [0, 100]
        If the VAF of a given variant for a given barcode
        is less than the given value, and the call is
        '1' i.e. HET, then it is converted to a no call.
    max_vaf : float [0, 100]
        If the VAF of a given variant for a given barcode
        is greater than the given value, and the call is
        '1' i.e. HET, then it is converted to a no call.
    min_prct_cells : float [0, 100]
        The minimum percent of total cells in which the variant
        should be called as '0', '1', or '2' after the
        filters are applied.
    min_mut_prct_cells : float [0, 100]
        The minimum percent of the total cells in which the
        variant should be mutated, i.e., either '1' or '2'
        after the filters are applied.
    min_std : float [0, 100]
        The standard deviation of the VAF across the cells
        of the variants should be greater than
        this value.

    Returns
    -------
    numpy.ndarray
    """

    gt = self.layers[NGT]
    dp = self.layers[DP]
    gq = self.layers[GQ]
    vaf = self.layers[AF]

    dp_keep = dp >= min_dp
    gq_keep = gq >= min_gq
    min_vaf_keep = ~np.logical_and(vaf < min_vaf, gt == 1)
    max_vaf_keep = ~np.logical_and(vaf > max_vaf, gt == 1)
    gt = (gt - 3) * dp_keep * gq_keep * min_vaf_keep * max_vaf_keep + 3  # workaround to apply filter in one line

    self.add_layer(NGT_FILTERED, gt)

    num_cells = len(self.barcodes())
    min_cells_filter = np.isin(gt, [0, 1, 2]).sum(axis=0) > num_cells * min_prct_cells / 100
    min_cells_mut_filter = np.isin(gt, [1, 2]).sum(axis=0) > num_cells * min_mut_prct_cells / 100
    good_variants = min_cells_mut_filter * min_cells_filter

    final_filter = (vaf.std(axis=0) >= min_std) * good_variants

    return self.col_attrs[ID][final_filter].astype(str)