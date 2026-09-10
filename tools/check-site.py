#!/usr/bin/env python3
"""Check generated local links/anchors and reject filesystem-only asset URLs."""
import json
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

root = Path(sys.argv[1] if len(sys.argv) > 1 else '_site').resolve()
class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.ids, self.links = set(), []
        self.feed(path.read_text(encoding='utf-8'))
    def handle_starttag(self, tag, attributes):
        a = dict(attributes)
        if 'id' in a:
            self.ids.add(a['id'])
        for key in ('href', 'src'):
            if a.get(key):
                self.links.append(a[key])

pages = {p.resolve(): Page(p) for p in root.rglob('*.html')}
failures, checked = [], 0
for path, page in pages.items():
    for link in page.links:
        url = urlsplit(link)
        if url.scheme == 'file':
            failures.append(f'{path.relative_to(root)}: filesystem URL {link}')
            continue
        if url.scheme or url.netloc or link.startswith('//'):
            continue
        target = (root / unquote(url.path).lstrip('/') if url.path.startswith('/')
                  else path.parent / unquote(url.path)).resolve() if url.path else path
        if target.is_dir():
            target /= 'index.html'
        checked += 1
        if not target.is_file():
            failures.append(f'{path.relative_to(root)}: missing {link}')
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            failures.append(f'{path.relative_to(root)}: missing anchor {link}')
print(json.dumps({'status': 'FAIL' if failures or not pages else 'PASS',
                  'html_pages': len(pages), 'local_links_checked': checked,
                  'failures': failures}, indent=2))
sys.exit(bool(failures or not pages))
