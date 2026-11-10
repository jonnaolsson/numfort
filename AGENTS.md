# Project instructions

## CI permission

- Do not create Docker images without explicit user approval of the proposed
  image and build commands, including temporary or debugging images.
- Obtain explicit user approval for the proposed commands and job scope before
  invoking the CI framework, including `--list` and container/image operations.
- Run the full matrix only when approval explicitly covers it.
- If approval has not been given, report validation as not run.

## Build and test validation

- For isolated debugging that does not require an unavailable gfortran version,
  use a temporary directory on the host instead of creating a Docker image.
- For CI validation, use the container-based framework in `ci/ubuntu26.04/`; follow
  [its README](ci/ubuntu26.04/README.md) for setup and commands.
- Use `run-tests.sh` for a single compiler/configuration or `run-matrix.sh` for
  selected named jobs.
- Preserve `build-ci` logs and matrix summaries. Check CTest results and
  `Testing/Temporary/LastTest.log`; executable exit status determines failure.

## Unit tests

- Unit tests are implemented in `tests/` and are run by CTest.
- Call each test subroutine directly from the main `program` body. Do not use
  intermediate wrapper subroutines such as `test_all()`.
- When writing new tests, split distinct test families into separate
  subroutines. For example, if testing `unique()`, combine tests for `unique()`
  into a single subroutine `test_unique()`, but do not combine tests for
  `unique()` and `sort()`.
- If a test family contains more than 10 test cases, split it into separate
  subroutines called directly from the main `program` body.
- Do not write tests that replicate the implementation verbatim.
- Create tests for known corner cases or closed-form solutions where the
  expected return value can be computed independently of the implementation.
- For deterministic algorithms with known reference implementations, verify
  correctness against hard-coded results from NumPy or SciPy.
- Use `error stop "<message>"` to terminate a test as soon as a failure is
  detected. Use static messages for error stops so that no boilerplate dynamic
  string creation is needed.
- Add line comments documenting each test case and data fixture.
- Add a header comment for each test subroutine with a one-sentence summary of
  its scope, followed by brief bullet points for each test case covered. Use
  concise descriptions and avoid repetitive lead-ins (such as "Test with...").
  Do not add separate bullets for iterations or parameter sweeps within a test
  case:
  ```fortran
  subroutine test_unique()
      !* Test unique() deduplication and array bounds.
      !   - Integer array with duplicates.
      !   - Floating-point array with duplicates.
      !   - Array sizes of 0, 1, 5, and 10.
  end subroutine
  ```

## Fortran style

- Follow the global Fortran style guide for modern code in numfort.
- For legacy external code in `src/external`, follow the style of the original
  code, but do not introduce obsolete constructs when adding new code to those
  files.
- Add brief line comments explaining non-trivial mathematical expressions,
  formulas, or derivations.

## Git

- Do not create commits without explicit user approval.
- Format commit subjects as `<prefix>: <summary>`.
- Use the modified numfort component as `<prefix>` (e.g., `linalg`, `optimize`),
  or an area like `cmake`, `tests`, `ci`, or `docs` if not component-specific.
- Include a bullet list of brief additional details in the commit body.

## Language and structure

- Keep explanations in Markdown and code comments brief and easy to read.
- Prefer short bullet points over long paragraphs.
- Organize Markdown documents into sections with descriptive headings.

## Main README audience

- Write the root `README.md` for numfort users, including those building from source.
- Keep internal details and instructions for developing numfort out of this README.

## User-facing documentation

- Describe only the current version.
- When refactoring or updating the CMake build system, remove outdated information
  and do not retain or add references to the previous version.
