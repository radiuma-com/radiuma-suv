"""Prevent distribution/import/API spelling drift across packaging and docs."""
from pathlib import Path
import ast
import re
import tomllib

ROOT = Path(__file__).resolve().parents[1]

def test_distribution_cli_and_exports():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert config["project"]["name"] == "radiuma-suv"
    assert config["project"]["scripts"]["radiuma-suv"] == "radiuma_suv.cli:main"
    initializer = ast.parse((ROOT / "src/radiuma_suv/__init__.py").read_text(encoding="utf-8"))
    assign = next(n for n in initializer.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "__all__" for t in n.targets))
    assert set(ast.literal_eval(assign.value)) == {"radiuma_suv", "radiuma_suv_from_input", "RadiumaSUVLoader"}
    core = ast.parse((ROOT / "src/radiuma_suv/core.py").read_text(encoding="utf-8"))
    declared = {n.name for n in core.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    assert {"radiuma_suv", "radiuma_suv_from_input", "RadiumaSUVLoader"} <= declared

def test_no_invalid_hyphenated_python_identifiers_in_docs():
    docs = [ROOT / "README.md", *(ROOT / "docs").glob("*.md"), *(ROOT / "examples").glob("*.py")]
    errors = []
    for doc in docs:
        content = doc.read_text(encoding="utf-8")
        if re.search(r"from\s+radiuma-suv\s+import|radiuma-suv_from_input|radiuma-suv\s*\(", content):
            errors.append(str(doc.relative_to(ROOT)))
    assert errors == []
