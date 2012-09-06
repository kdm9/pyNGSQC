def _percentile_from_counts(count_list, percentile):
    """Returns the median value from a list whose values are counts of the
    index i
    """
    # The index of the median, or of the lower of the two values if sum is even
    halfway_index = int(round(float(sum(count_list) * percentile))) - 1
    current_index = 0
    pos_within_value = 0  # Governs when we skip to the next median value
    median = 0
    while current_index <= halfway_index:
        # If we have not iterated through all counts of this value of median
        if  pos_within_value < count_list[median]:
            pos_within_value += 1
            current_index += 1
        else:
            # Move to next median
            median += 1
            pos_within_value = 0

    # At this point, median == lower of two values if the number of counts is
    # even, so we need to average the current value of median with the
    # following value of median (which may be the same number) to get the
    # median. If the number of counts is odd, then median == true median
    if sum(count_list) % 2 == 0:
        lower_score = median
        if count_list[median] < pos_within_value + 1:
            upper_score = median + 1
        else:
            upper_score = median
        median = float(lower_score + upper_score) / 2.0
    return float(median)

l = [16511,1651,51,516,516,156,51,51,6516,516,516,15,1,61,65,16,51,651,615,65,16,51,6151,65,165,16,5161,651,65,16,516,51,651,65,16,516,51,651,51615,1651,651,6516,516,51]
for p in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
    print _percentile_from_counts(l, p)
