"""Where the page exporters put their frames and the builder its pages."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
DATA = ROOT / "analysis" / ".pages" / "data"
PAGES = ROOT / "analysis" / ".pages"
TEMPLATES = Path(__file__).resolve().parent / "templates"


def out(name):
    DATA.mkdir(parents=True, exist_ok=True)
    return DATA / f"{name}.json"
