def sanitize_3p(value: str):
    if value.startswith("^"):
        raise AssertionError("A 3p adapter must not start with ^")

    return value

def sanitize_5p(value: str):
    if value.endswith("$"):
        raise AssertionError("A 5p adapter must not end with $")

    return value