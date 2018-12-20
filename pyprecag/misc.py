def roundDownToPixel(num, divisor):
    """ Round a number DOWN to the nearest divisor (ie pixel size)
        eg. 606131.297544 -> divisor=2.5 -> 606130.0
    This ensures a floating number is a multiple of the divisor (ie pixel size) by rounding it down. This is useful for
    getting the min extent of a vector dataset
     can be used to
    determine the minimum X or Y coordinate of a pixel edge across multiple datasets to ensure the alignment of pixels
    when converting vectors to rasters.

    Args:
        num (float): A number which could represent a
        divisor (float):

    Returns:
        float:

    Examples:
        >>> roundDownToPixel(606131.297544,2.5)
        606130.0
    """
    return num - (num % divisor)


def roundUpToPixel(num, divisor):
    """Round a number UP to the nearest divisor
       eg 607426.0264498683 -> divisor=2.5 -> 607427.5
    Args:
        num (float):
        divisor (float):
    Returns:
        float:

    Examples:
        >>> roundUpToPixel(607426.0264498683 , 2.5)
        607427.5

    """
    num = num + divisor
    return num - (num % divisor)
