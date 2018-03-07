
def remove_softclipping(read, cigar):
    pos = 0
    buf = ""
    new_read = ""
    new_cigar = ""

    for char in cigar:
        if char in "0123456789":
            buf += char
            continue

        length = int(buf)

        if char != "S":
            if char != "N":
                new_read += read[pos:pos+length]
            new_cigar += buf + char

        if char != "N":
            pos += length
        buf = ""

    return new_read, new_cigar