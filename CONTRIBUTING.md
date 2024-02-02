Contributing
============

Contributions to the code are highly welcome.
Please consider the
[Treelm Coding Guidelines](https://geb.inf.tu-dresden.de/doxy/treelm/page/codingGuidelines.html)
for the contributed code.

Please note the complication by the split into
submodules.
The source code of the solver resides in the
submodule musubi-source, which is linked to here
in the `mus` subdiectory.

To contribute code, please open a new branch in
both, the repository in `mus` and here in the
parent repository, i.e.:

```
cd mus
git checkout -b <yourbranchname>
cd ..
git checkout -b <yourbranchname>
```

When you have changes on the branch and want to
feed it back, open a pull request for both branches.
The pull request in the `mus` submodule will have
to be resolved first and merged back.
With that in place we then proceed with the merge
request in the parent repository.
The merged change should always only refer to commits
on the main branch of the submodules, and thus, the
merges here that update references to the submodules
are always squashed.
For the pull request in the parent directory, use the
`?template=sub_pr.md` query parameter in the URL for
the pull request, that is:

```
https://github.com/apes-suite/musubi/compare/<branch>?template=sub_pr.md
```
