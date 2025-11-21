import re
import unicodedata

def safe_filename(name, replacement="_", max_length=255):
    # Normalize unicode (é → e, etc.)
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()

    # Replace hyphens first
    name = name.replace("-", "_")

    # Replace any remaining invalid chars
    name = re.sub(r"[^A-Za-z0-9._-]", replacement, name)

    # Remove repeated replacement characters
    name = re.sub(r"_+", "_", name)

    # Trim to safe max filename length
    return name