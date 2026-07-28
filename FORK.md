# Fork of Oligo Designer Toolsuite

This file more or less explains how this fork works. All emdashes are human-made – even this one.

## Git Setup

The primary development branch on this fork is `bachelor-thesis`. Direct commits are allowed for small changes – anything that might be relevant upstream needs a PR.
`bachelor-thesis` is branched off `dev`, which is the precursor to `main` (at least upstream). We will mostly ignore `main` here.

We will sync any changes in upstream's `dev` into our `dev`. They will then get merged into our `bachelor-thesis`. Or `bachelor-thesis` will be rebased onto `dev`, not decided yet.

## Contributing upstream

The focus of this fork is to be the codebase supporting a bachelor thesis, not to directly contribute changes and improvements upstream. This allows us more creative freedom in how we manage the project (especially i.r.t. tooling and dependencies).

Nonetheless, we would love to see the improvements made here go upstream, so we don't want to make it needlessly difficult to integrate them after-the-fact. For example, we still use PRs for isolated changes and try to keep the full commit history for retroactive integration. We might even keep all merged branches if that doesn't become a mess.
