---
name: branch-update
description: Use proactively when reviewing code to update the source branch to be current with its target branch while maintaining proper git history (avoiding foxtrot merges).
tools: Bash(git checkout *), Bash(git merge *), Bash(git branch *), Bash(git status), TodoWrite
---

**When to Use:**

- When your feature branch needs to be updated with the latest changes from the
  target branch
- When preparing a branch for review and it's behind the target branch
- When you want to ensure the target branch appears as the left parent in merge
  history

**Prerequisites:**

- Working directory must be clean (all changes committed and pushed to origin)
- Source branch should exist and be pushed to origin
- Target branch should be specified or determinable

**Process:**

1. **Validate Prerequisites:**

   - Check git status is clean
   - Verify source branch exists
   - Confirm target branch exists
   - Ensure current changes are pushed to origin

2. **Create Temporary Branch:**

   - `git checkout -b <source-branch-name>-temp <target-branch>`
   - This creates a temporary branch starting from the target branch

3. **Merge Source Branch:**

   - `git merge <source-branch-name>`
   - Merge the source branch into the temporary branch

4. **Handle Conflicts:**

   - If conflicts occur, provide guidance to user
   - Wait for user to resolve conflicts and commit
   - Continue once conflicts are resolved

5. **Update Source Branch:**

   - `git branch -f <source-branch-name>` (force update source branch to current
     HEAD)
   - `git checkout <source-branch-name>` (switch back to source branch)
   - `git branch -d <source-branch-name>-temp` (delete temporary branch)

6. **Verify Results:**
   - Confirm the branch update was successful
   - Show the updated branch status

**Benefits:**

- Maintains clean git history with target branch as left parent
- Avoids foxtrot merges that can complicate history
- Preserves original commit structure while bringing in target changes
- Safe process that can be aborted if issues arise

**Error Handling:**

- If working directory is not clean, abort with instructions
- If branches don't exist, provide clear error messages
- If merge conflicts occur, guide user through resolution
- If any step fails, provide recovery instructions

**Example Usage:**

```
# Starting state: feature branch behind main
# Goal: Update feature branch with latest main changes
# Result: Clean history with main as left parent
```

**Integration:** This agent can be used by review-targeting commands instead of
aborting when branches are out of sync, providing a seamless workflow for
keeping feature branches current.
