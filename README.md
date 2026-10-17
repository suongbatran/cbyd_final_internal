This repository contains ...

---

## Browse the routes: `dashboard.html`

All routes can be explored in an interactive viewer that is **one self-contained HTML file**.
It needs no installation, no server and no internet connection.

1. Download `[dashboard.html](dashboard.html)`.
2. Open it in any browser.

**What you can do in the viewer**

- **Search** routes by compound name, DOI, author or year.
- **See each route as a synthesis tree.** Starting materials are 🟢, intermediates 🔵,
targets 🔴 and reaction steps 🟨. Each structure is drawn in the same orientation as the
product it becomes, so you can follow the scaffold from step to step.
- **Click a reaction** to see the atom-mapped reaction, the original description, the
reagents and conditions (hover a reagent thumbnail to enlarge it) and the mapped SMILES.



## Found a mistake? Tell us

In the viewer, click **Report a mistake** (or open the report box under any
reaction or molecule). The report is pre-filled with the route, the step and its SMILES.
Then choose one of these:


| Option                 | What happens                                                                                                                         |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| **Email**              | Opens your mail app with a pre-filled message to $\color{red}{\textsf{xxx}}$                                                         |
| **Open GitHub issue**  | Opens a pre-filled [issue](../../issues/new) on this repository                                                                      |
| **Suggest a fix (PR)** | Opens the route's JSON file in GitHub's editor. Make the change, click **Propose changes**, and GitHub opens a pull request for you. |


You can also edit any file in `[data/routes/](data/routes)` directly and open a pull request.

---



## Data format

There is one JSON file per route in `[data/routes/](data/routes)`, named `<Compound><Year><Author>_processed.json`.

```jsonc
{
  "background": {
    "name": ...,                     // target compound
    "author": ...,
    "year": ...,
    "drug name": ...,                // set for marketed-drug routes
    "manufacturer": ...,
    "citation": ...,
    "doi": ...
  },
  "synthesis": [                     // one entry per reaction step
    {
      "reaction_id": ...,
      "raw_description": ...,        // conditions as reported
      "reaction_smiles": ...,
      "full_reaction_smiles": ...,
      "full_reaction_smiles_mapped": ...,
      "mapping_confidence": ...,
      "parsed_reaction": {
        "GENERAL_REAGENTS": [...],  "GENERAL_REAGENTS_DEFINITIONS": [...],
        "CATALYSTS": [...],         "CATALYSTS_DEFINITIONS": [...],
        "LIGANDS": [...],           "LIGANDS_DEFINITIONS": [...],
        "SOLVENTS": [...],          "SOLVENTS_DEFINITIONS": [...],
        "MODIFIERS": [...],
        "TEMPERATURE": [...],
        "YIELD": [...]
      }
    }
  ]
}
```

- `reaction_smiles`: the main transformation (reactants `>>` product).
- `full_reaction_smiles`: the same reaction with reagents added from `parsed_reaction`.
- `full_reaction_smiles_mapped`: the atom-mapped reaction, with `mapping_confidence` from the mapper (RXNMapper).
- Each `*_DEFINITIONS` list holds the SMILES for the names in the list with the same prefix, in the same order.



## Repository layout


| Path             | Contents                                                                                      |
| ---------------- | --------------------------------------------------------------------------------------------- |
| `dashboard.html` | The standalone viewer: all routes, all code, works offline                                    |
| `data/routes/`   | The dataset, one JSON file per route                                                          |
| `dashboard/`     | Source code of the viewer. `build_static.py` regenerates `dashboard.html` from `data/routes/` |




## Citation

If you use this dataset, please cite:

```bibtex
@article{xxx,
  title   = {xxx},
  author  = {xxx},
  journal = {xxx},
  year    = {xxx},
  doi     = {xxx}
}
```

