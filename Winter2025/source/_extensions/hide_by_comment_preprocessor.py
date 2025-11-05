from nbconvert.preprocessors import Preprocessor
import sys

class HideByCommentPreprocessor(Preprocessor):
    def preprocess_cell(self, cell, resources, index):
        src = cell.get("source") or ""
        if not src:
            return cell, resources
        first = src.split("\n", 1)[0].strip()
        # match '# hide' at the first line (case-insensitive, spaces allowed)
        if first.startswith("#") and "hide" in first.lower():
            print(f"[HBCP] hiding cell #{index}", file=sys.stderr)
            # keep only the comment line; to remove the whole cell, return (None, resources)
            cell["source"] = f"{first}: the code in this cell is hidden by the author"
        return cell, resources
