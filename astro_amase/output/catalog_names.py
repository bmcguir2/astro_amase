"""Keep scientific catalog labels intact while writing portable output filenames."""
from collections import Counter
from hashlib import sha256
import unicodedata


def catalog_file_stems(labels):
    """Map display labels to distinct stems, including on case-insensitive filesystems."""
    def escape(label):
        return "".join(
            f"%{ord(char):02X}" if char in '%~/\\<>:"|?*' or ord(char) < 32 or ord(char) == 127 else char
            for char in label
        )

    def filesystem_key(name):
        return unicodedata.normalize("NFD", name).casefold()

    stems = {label: escape(label) for label in labels}
    counts = Counter(filesystem_key(stem) for stem in stems.values())
    for label, stem in stems.items():
        key = filesystem_key(stem)
        if counts[key] > 1 or key == "total" or len(stem.encode("utf-8")) > 205 or not stem:
            stems[label] = "catalog~" + sha256(label.encode("utf-8")).hexdigest()
    return stems
