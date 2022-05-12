"""
Original file is located at
    https://colab.research.google.com/drive/17K7e9YIkcJ8dcNsSRUEXw7KxU8zPtdO_
"""
import math
import numpy
import scipy.special as spc


# -----------Random Generator-----------
def lfsr(seed, taps, nbits):
    """
    seed: binary value
    taps: (tuple) every value in the tuple represents a bit in the seed
    nbits: number of bits
    """
    sr = seed
    while 1:
        xor = 1
        for t in taps:
            if (sr & (1 << (t - 1))) != 0:
                xor ^= 1
        sr = (xor << nbits - 1) + (sr >> 1)
        yield xor, sr
        if sr == seed:
            break


# --------------Testers------------------
def monobit(bin_data: str):
    """
    Note that this description is taken from the NIST documentation [1]
    [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf

    The focus of this test is the proportion of zeros and ones for the entire sequence. The purpose of this test is
    to determine whether the number of ones and zeros in a sequence are approximately the same as would be expected
    for a truly random sequence. This test assesses the closeness of the fraction of ones to 1/2, that is the number
    of ones and zeros ina  sequence should be about the same. All subsequent tests depend on this test.

    :param bin_data: a binary string
    :return: the p-value from the test
    """
    count = 0
    # If the char is 0 minus 1, else add 1
    for char in bin_data:
        if char == "0":
            count -= 1
        else:
            count += 1
    # Calculate the p value
    sobs = count / math.sqrt(len(bin_data))
    p_val = spc.erfc(math.fabs(sobs) / math.sqrt(2))
    return p_val


def block_frequency(bin_data: str, block_size=128):
    """
    Note that this description is taken from the NIST documentation [1]
    [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf
    The focus of this tests is the proportion of ones within M-bit blocks. The purpose of this tests is to determine
    whether the frequency of ones in an M-bit block is approximately M/2, as would be expected under an assumption
    of randomness. For block size M=1, this test degenerates to the monobit frequency test.
    :param bin_data: a binary string
    :return: the p-value from the test
    :param block_size: the size of the blocks that the binary sequence is partitioned into
    """
    # Work out the number of blocks, discard the remainder
    num_blocks = math.floor(len(bin_data) / block_size)
    block_start, block_end = 0, block_size
    # Keep track of the proportion of ones per block
    proportion_sum = 0.0
    for i in range(num_blocks):
        # Slice the binary string into a block
        block_data = bin_data[block_start:block_end]
        # Keep track of the number of ones
        ones_count = 0
        for char in block_data:
            if char == "1":
                ones_count += 1
        pi = ones_count / block_size
        proportion_sum += pow(pi - 0.5, 2.0)
        # Update the slice locations
        block_start += block_size
        block_end += block_size
    # Calculate the p-value
    chi_squared = 4.0 * block_size * proportion_sum
    p_val = spc.gammaincc(num_blocks / 2, chi_squared / 2)
    return p_val


def serial(bin_data, pattern_length=16, method="first"):
    """
    Note that this description is taken from the NIST documentation [1]
    [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf
    The focus of this test is the frequency of all possible overlapping m-bit patterns across the entire
    sequence. The purpose of this test is to determine whether the number of occurrences of the 2m m-bit
    overlapping patterns is approximately the same as would be expected for a random sequence. Random
    sequences have uniformity; that is, every m-bit pattern has the same chance of appearing as every other
    m-bit pattern. Note that for m = 1, the Serial test is equivalent to the Frequency test of Section 2.1.
    :param bin_data: a binary string
    :param pattern_length: the length of the pattern (m)
    :return: the P value
    """
    n = len(bin_data)
    # Add first m-1 bits to the end
    bin_data += bin_data[: pattern_length - 1 :]

    # Get max length one patterns for m, m-1, m-2
    max_pattern = ""
    for i in range(pattern_length + 1):
        max_pattern += "1"

    # Keep track of each pattern's frequency (how often it appears)
    vobs_one = numpy.zeros(int(max_pattern[0:pattern_length:], 2) + 1)
    vobs_two = numpy.zeros(int(max_pattern[0 : pattern_length - 1 :], 2) + 1)
    vobs_thr = numpy.zeros(int(max_pattern[0 : pattern_length - 2 :], 2) + 1)

    for i in range(n):
        # Work out what pattern is observed
        vobs_one[int(bin_data[i : i + pattern_length :], 2)] += 1
        vobs_two[int(bin_data[i : i + pattern_length - 1 :], 2)] += 1
        vobs_thr[int(bin_data[i : i + pattern_length - 2 :], 2)] += 1

    vobs = [vobs_one, vobs_two, vobs_thr]
    sums = numpy.zeros(3)
    for i in range(3):
        for j in range(len(vobs[i])):
            sums[i] += pow(vobs[i][j], 2)
        sums[i] = (sums[i] * pow(2, pattern_length - i) / n) - n

    # Calculate the test statistics and p values
    del1 = sums[0] - sums[1]
    del2 = sums[0] - 2.0 * sums[1] + sums[2]
    p_val_one = spc.gammaincc(pow(2, pattern_length - 1) / 2, del1 / 2.0)
    p_val_two = spc.gammaincc(pow(2, pattern_length - 2) / 2, del2 / 2.0)

    # For checking the outputs
    if method == "first":
        return p_val_one
    else:
        # I am not sure if this is correct, but it makes sense to me.
        return min(p_val_one, p_val_two)


# ---------------executing tests--------------------
nbits = 8
l = []
for xor, sr in lfsr2(0b10011011, (8, 7, 6, 1), nbits):
    l.append(xor)
# The tests requires a string, so is necesary to convert the int list to string.
l_string = "".join(str(e) for e in l)

print("\nMonobit:", monobit(l_string))
print("Frequency Block:", block_frequency(l_string))
print("Serial:", serial(l_string))
print("\n_______________________\n")
print(
    "If the result of the test is major than 0.1 then the generator is random. Else, it's not random."
)
print("_______________________")
