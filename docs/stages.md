# Preprocessing Stages

## Overview

Preprocessing is divided into independent stages.  
Each stage writes to its own directory and saves a tagged `.set` file.

---

## Stage Order
1. Filter
2. Notch
3. Resample
4. Re-reference
5. INITREJ
6. ICA
7. ICLabel
8. Epoch
9. Baseline

---

## Stage Philosophy
Stages are:
- Explicit
- Isolated
- Reproducible
- Inspectable

No stage overwrites previous outputs.

---

## Resume Behaviour
Before executing a stage:
- The pipeline checks if the tagged `.set` file exists.
- If it exists → it is loaded
- If not → it is computed

This allows:
- Recovery after crashes
- Skipping completed steps
- Safe long-running preprocessing

---

## Manual Stages
Two stages involve manual decisions:
- INITREJ (channel interpolation)
- IC rejection

All manual decisions are logged.  
No automatic removal occurs without explicit user confirmation.