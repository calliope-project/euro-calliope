def unit(value, unit, parenthesis=True):
    """Formats the numeric value of a unit into a string in a consistent way."""
    formatted = f"{value:,g} {unit}"
    if parenthesis:
        formatted = f"({formatted})"
    return formatted
