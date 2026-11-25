"""Audit build paths against Git's case-sensitive tree, even on a Mac filesystem."""
import ast
from collections import defaultdict
from pathlib import Path, PurePosixPath
import re
import subprocess
import unittest

ROOT=Path(__file__).resolve().parents[2]


class BuildPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tracked=set(subprocess.check_output(['git','-C',str(ROOT),'ls-files'],text=True).splitlines())

    def test_scons_source_paths_match_tracked_spelling(self):
        # Path.exists() would miss this regression on case-insensitive mounts.
        tree=ast.parse((ROOT/'source/SConstruct').read_text())
        sources=set()
        for node in ast.walk(tree):
            if isinstance(node,ast.Constant) and isinstance(node.value,str):
                for token in node.value.split():
                    if re.fullmatch(r'[\w./-]+\.(?:cc|cpp|c)',token):
                        sources.add('source/'+str(PurePosixPath(token)))
        self.assertGreater(len(sources),80)
        missing=sorted(sources-self.tracked)
        self.assertEqual(missing,[],f'SCons sources not present with exact Git spelling: {missing}')

    def test_project_header_include_case(self):
        headers=defaultdict(set)
        for name in self.tracked:
            if name.endswith(('.h','.hh','.hpp')):
                base=PurePosixPath(name).name;headers[base.lower()].add(base)
        mismatches=[]
        for name in sorted(self.tracked):
            if not name.startswith('source/') or not name.endswith(('.cc','.h')):continue
            text=(ROOT/name).read_text(errors='replace')
            for include in re.findall(r'^\s*#\s*include\s*"([^"]+)"',text,re.M):
                base=PurePosixPath(include).name;spellings=headers.get(base.lower(),set())
                if spellings and base not in spellings:mismatches.append((name,include,sorted(spellings)))
        self.assertEqual(mismatches,[])


if __name__=='__main__':unittest.main()
