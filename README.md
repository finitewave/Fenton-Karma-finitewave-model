## Finitewave model template (replace with the model name)

Add model description here and fill the sections below.

This model implementation can be used separately from the Finitewave, allowing for standalone simulations and testing of the model dynamics without the need for the entire framework.

### Developer TODO Checklist (template repository only)
- [ ] **Change model entry point in `pyproject.toml`**  
  Update the entry point in the `pyproject.toml` file to reflect the new model's identifier. This ensures that the model can be correctly referenced and utilized within the project.
  
  In pyproject.toml, replace "model_template" with the actual model id. It must match the name of the directory where ops.py is located.
  ```toml
  [project.entry-points."finitewave.models"] 
  model_template = "finitewave_models.model_template"
  ```

- [ ] **Implement `ops.py` model equations**  
  Implement the model equations in the `ops.py` module. This module is the single source of truth for the model equations. Provide pure Python functions with scalar inputs/outputs (no NumPy arrays, no classes, no globals). Do NOT add numba/jax/torch here — the Finitewave runtime will wrap these functions for you. Stimulus and time integration are handled outside of the model. Here you only return time derivatives.

- [ ] **Create 0D implementation for tests and examples**  
  Design a zero-dimensional (0D) implementation of the model that can be used for testing and demonstration purposes. This should simplify the model to its core functionality, making it easier to validate and showcase.

- [ ] **Provide at least one workable example for the model**  
  Create a practical example that demonstrates how to use the model: model initialization, model stimulation and visualization of the model's AP.

- [ ] **Implement model tests**  
  Write unit tests to verify the correctness of the model's implementation: check model attributes and estimate AP parameters (min/max amplitude, duration) - this ensures that the model is excitable and generates the expected AP. 

### Reference
Paper, Authors, DOI.

### How to use (quickstart)
```bash
python -m example/model_example
```

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── model_template/                  # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       # fill with the model equations (pure functions)
├── implementation/                  # 0D model implementation
│   ├── __init__.py
│   └── model_0d.py
├── example/
│   └── model_example.py             # minimal script to run a short trace
├── tests/
│   └── test.py                      # smoke test; extend with reproducibility checks
├── .gitignore
├── LICENSE                          # MIT
├── pyproject.toml                   # placeholders to replace
└── README.md                        # this file
```

### Variables
Model state variables: description, units and ranges (optional)
- `u` — ...

### Parameters
Parameters and their defualt values
- `par` - ...

