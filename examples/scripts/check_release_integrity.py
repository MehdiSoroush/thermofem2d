"""Validate local release documentation references and version metadata.

This prepublication check addresses structural completeness, not numerical accuracy.
Run from the project root with ``python examples/scripts/check_release_integrity.py``.
"""
from pathlib import Path
import json
import re
import tomllib


ROOT = Path(__file__).resolve().parents[2]
MARKDOWN_LINK = re.compile(r'(?<!\!)\[[^\]\n]+\]\(([^)]+)\)')
IMAGE_LINK = re.compile(r'!\[[^\]\n]*\]\(([^)]+)\)')


def local_links(document):
    text = document.read_text(encoding='utf-8')
    for target in [*MARKDOWN_LINK.findall(text), *IMAGE_LINK.findall(text)]:
        target = target.split('#', 1)[0].strip()
        if not target or target.startswith(('https://', 'http://', 'mailto:')):
            continue
        yield target


def check_release():
    issues = []
    docs = [ROOT / 'README.md', *sorted((ROOT / 'docs').glob('*.md'))]
    for doc in docs:
        for target in local_links(doc):
            if not (doc.parent / target).exists():
                issues.append(f'{doc.relative_to(ROOT)}: missing {target}')
    pyproject = tomllib.loads((ROOT / 'pyproject.toml').read_text(encoding='utf-8'))
    version = pyproject['project']['version']
    release_notes = ROOT / 'docs' / f'RELEASE_NOTES_{version}.md'
    if not release_notes.is_file():
        issues.append(f'missing release notes: {release_notes.relative_to(ROOT)}')
    changelog = (ROOT / 'CHANGELOG.md').read_text(encoding='utf-8')
    if not re.search(rf'^## v?{re.escape(version)}(?:\\s|$)', changelog, flags=re.MULTILINE):
        issues.append(f'CHANGELOG.md: missing version heading {version}')
    package = (ROOT / 'src/thermofem2d/__init__.py').read_text(encoding='utf-8')
    citation = (ROOT / 'CITATION.cff').read_text(encoding='utf-8')
    readme = (ROOT / 'README.md').read_text(encoding='utf-8')
    for file,pattern in [
        ('__init__.py', f'__version__ = "{version}"'),
        ('CITATION.cff', f'version: {version}'),
        ('README.md', f'v{version}'),
    ]:
        source = {'__init__.py':package,'CITATION.cff':citation,'README.md':readme}[file]
        if pattern not in source:
            issues.append(f'{file}: inconsistent release version')
    for path in ROOT.glob('PRIVATE_*'):
        issues.append(f'private transfer material present: {path.name}')
    for path in [*sorted((ROOT/'results/verification').glob('v07_*.json'))]:
        data=json.loads(path.read_text(encoding='utf-8'))
        if any(w in json.dumps(data).lower() for w in
               ('publication_status": "hold','not a public release')):
            issues.append(f'{path.relative_to(ROOT)}: stale release metadata')
    if issues:
        raise SystemExit('Release integrity check failed:\n'+'\n'.join(issues))
    print(f'Release metadata and {len(docs)} Markdown documents checked: OK (v{version})')


if __name__ == '__main__':
    check_release()
