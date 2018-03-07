import os
import subprocess

from .exceptions import QurallkyException
from .conf_sanitation import adapter, general

_sanitation = {
    "adapters": {
        "3p": adapter.sanitize_3p,
        "5p": adapter.sanitize_5p,
    },
    "general": {
        "input": (general.directory, {}),
        "output": (general.directory, {}),
        "threads": int,
    },
    "cutadapt": {
        "minReadLength": int,
    },
    "modcount": {
        "minReadsPerGene": int,
    }
}

class Configuration:
    """
    Reads in configuration.
    """
    config = {}
    samples = []
    force_overwrite = False

    def __init__(self, filename: str):
        self.config = {}
        self.samples = []

        self.base_path = os.path.dirname(filename)
        self.read_config(filename)
        self.read_samples()

    def read_config(self, filename: str):
        if not os.path.exists(filename):
            raise ConfigurationException("Configuration file {} does not exist.")
        if not os.path.isfile(filename):
            raise ConfigurationException("Configuration file {} is not a valid file.")

        with open(filename, "r") as fh:
            i = 0
            section = None
            for line in fh:
                i+=1
                line = line.strip()

                if len(line) == 0:
                    continue

                # # at the beginning of a line marks a comment, we skip those.
                if line.startswith("#"):
                    continue

                if line.startswith("["):
                    if line.endswith("]") is False:
                        raise ConfigurationException(
                            "A section marker has been opened but not closed (line {}).".format(i))

                    section = line[1:-1].strip()
                    continue

                if section is not None:
                    keyval = [x.strip() for x in line.split("=", 1)]

                    if len(keyval) == 2:
                        key, val = keyval
                    else:
                        key, val = keyval, None
                self.set_config(section, key, val)

    def set_config(self, section: str, key: str, val: str):
        if section in _sanitation and key in _sanitation[section]:
            try:
                sanitize = _sanitation[section][key]
                if isinstance(sanitize, tuple) or isinstance(sanitize, list):
                    val = sanitize[0](val, **{**sanitize[1], **{"base_dir": self.base_path}})
                else:
                    val = sanitize(val)
            except AssertionError as e:
                print(e)
                raise ConfigurationException(str(e))

        if section not in self.config:
            self.config[section] = {}
        self.config[section][key] = val

    def get(self, section: str, key: str):
        if section in self.config and key in self.config[section]:
            return self.config[section][key]
        else:
            raise ConfigurationException("Configuration.get: {}.{} has not been found".format(section, key))

    def read_samples(self):
        for file in os.listdir(self.get("general", "input")):
            if file.startswith("."):
                continue

            if file.endswith("fastq.gz") is False and file.endswith(".fastq") is False:
                continue

            self.samples.append(file)

    def extend_path(self, filename):
        return os.path.join(self.base_path, filename)

    def extend_input_path(self, filename):
        return os.path.join(self.get("general", "input"), filename)

    def extend_output_path(self, filename):
        return os.path.join(self.get("general", "output"), filename)

    def execute(self, *commands):
        for command in commands:
            subprocess.check_output(command, shell=True)
            #print(command)


class ConfigurationException(QurallkyException):
    pass

class DoContinue():
    pass