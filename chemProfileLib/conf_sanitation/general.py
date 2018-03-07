import os

def directory(val: str, base_dir: str, **kwargs) -> str:
    if val.startswith("/"):
        return val
    else:
        return os.path.join(base_dir, val)