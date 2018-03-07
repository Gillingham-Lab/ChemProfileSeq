# !/usr/bin/env python
import argh
import sys

# import commands
from chemProfileLib import QurallkyException
from chemProfileLib.cli import name, version, t, debug
from chemProfileLib.cli.toolchains import align, check, plotStops, logo, count, enrichment_table


def main():
    try:
        print(t.bold("{} version {}\n").format(name, version))

        if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ["-h", "help", "--help"]):
            print(
                "ChemProfileSeq toolchain for analysing chemically induced stops or mutations from NGS sequencing data.\n\n")

        argh.dispatch_commands(sorted([
            align, check,
            count, enrichment_table,
            plotStops, logo,
        ], key=lambda x: x.__name__))

    except KeyboardInterrupt:
        print("\nInterrupted; Data might be corrupted.")
    except QurallkyException as e:
        print(t.red("An error occured:"))
        print("\t", e)
    except Exception as e:
        print(e)
        if debug:
            raise e