"""Check shipped notebook, analysis source, and the offline instructions using only Python's standard library."""
import ast
from base64 import b64decode
from hashlib import sha256
from html.parser import HTMLParser
from io import BytesIO
import json
from pathlib import Path
from urllib.parse import urlsplit
import unittest
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.anchors = []
        self.downloads = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.anchors.append(attrs["id"])
        for name in ("href", "src"):
            if name in attrs:
                self.links.append((tag, attrs[name]))
                if tag == "a" and "download" in attrs:
                    self.downloads.append(attrs[name])


class LiteDistribution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        name = "downloads/AstroAmase_Lite_Integrated.ipynb"
        cls.notebook = json.loads((DOCS / name).read_bytes())
        cls.page = Page()
        cls.page.feed((DOCS / "index.html").read_text(encoding="utf-8"))
        cls.cells = {c["id"]: c for c in cls.notebook["cells"] if c["cell_type"] == "code"}
        cls.assets = {}
        wanted = {"ENGINE_RECORD", "ENGINE_BASE64", "ENGINE_SHA", "PUBLIC_RECIPES", "NATIVE_RECIPE"}
        tree = ast.parse("".join(cls.cells["setup-environment"]["source"]))
        for node in tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                name = node.targets[0].id
                if name in wanted:
                    cls.assets[name] = ast.literal_eval(node.value)

    def test_clean_executable_notebook(self):
        self.assertEqual(self.notebook["nbformat"], 4)
        self.assertEqual(len({c["id"] for c in self.notebook["cells"]}), len(self.notebook["cells"]))
        for name, cell in self.cells.items():
            with self.subTest(cell=name):
                self.assertIsNone(cell["execution_count"])
                self.assertEqual(cell["outputs"], [])
                ast.parse("".join(cell["source"]))

    def test_engine_is_the_checked_out_source(self):
        record = self.assets["ENGINE_RECORD"]
        engine = b64decode(self.assets["ENGINE_BASE64"], validate=True)
        self.assertEqual(sha256(engine).hexdigest(), self.assets["ENGINE_SHA"])
        self.assertEqual(self.assets["ENGINE_SHA"], record["sha256"])
        with ZipFile(BytesIO(engine)) as bundle:
            self.assertIsNone(bundle.testzip())
            self.assertCountEqual(bundle.namelist(), list(record["source_files"]) + ["licenses/molsim-LICENSE"])
            self.assertEqual(bundle.read("licenses/molsim-LICENSE"), (ROOT / "licenses/molsim-GPL-3.0.txt").read_bytes())
            for name, digest in record["source_files"].items():
                with self.subTest(file=name):
                    data = bundle.read(name)
                    self.assertEqual(data, (ROOT/name).read_bytes())
                    self.assertEqual(sha256(data).hexdigest(), digest)

    def test_recipe_pins(self):
        for spec in [*self.assets["PUBLIC_RECIPES"].values(), self.assets["NATIVE_RECIPE"]]:
            raw = b64decode(spec["blob"], validate=True)
            self.assertEqual(sha256(raw).hexdigest(), spec["sha256"])

    def test_offline_guide(self):
        self.assertEqual(len(set(self.page.anchors)), len(self.page.anchors))
        self.assertIn("downloads/AstroAmase_Lite_Integrated.ipynb", self.page.downloads)
        self.assertIn("assets/rotor-upload-example.zip", self.page.downloads)
        for tag, target in self.page.links:
            with self.subTest(link=target):
                parts = urlsplit(target)
                if parts.scheme:
                    self.assertEqual(tag, "a")
                    self.assertEqual(parts.scheme, "https")
                elif parts.path:
                    self.assertTrue((DOCS/parts.path).is_file())
                elif parts.fragment:
                    self.assertIn(parts.fragment, self.page.anchors)
        with ZipFile(DOCS/"assets/rotor-upload-example.zip") as z:
            self.assertIsNone(z.testzip())
            self.assertEqual({Path(name).name for name in z.namelist()},
                             {"spectrum.txt","metadata.json","transitions.csv","injection.json"})

if __name__ == "__main__":
    unittest.main()
