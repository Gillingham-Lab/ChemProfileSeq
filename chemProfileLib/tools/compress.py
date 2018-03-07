import subprocess

def run(config, infile):
    compress = "gzip -f {infile}".format(infile=infile)
    config.execute(compress)

    return infile + ".gz"