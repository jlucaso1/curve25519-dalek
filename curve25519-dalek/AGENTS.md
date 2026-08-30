# Instructions for an agent

## Before writing code

* If you haven't already read the `AGENTS.md` at the root of this project, please do so now.

## When writing code

* This crate does not have `rand` as a dev dependency. It does have `proptest`, though. Thus, when writing tests that use randomness, prefer `proptest` over some ad-hoc randomness method
