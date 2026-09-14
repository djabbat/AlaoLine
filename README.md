# AlaoLine

**N₂-atmosphere short-sprouting for bioavailable iron from wheat.**

AlaoLine is an open-source software package for the modelling, design and
optimisation of the Alao N₂-atmosphere short-sprouting process: a circular
bio-based process that converts wheat into a bioavailable-iron food
ingredient.

The package implements:
- a kinetic model of Fe²⁺ preservation and IP₆ (phytate) hydrolysis across
  the five processing stages;
- the 2⁽⁵⁻¹⁾ fractional-factorial experimental design (resolution V) used in
  WP3;
- response-surface (RSM) fitting and optimisation.

It is released open-source (MIT) as part of the CIRCBIO-07 AlaoLine proposal
(HORIZON-CL6-2026-01-CIRCBIO-07), coordinated by Georgia Longevity Alliance.

## Install

```bash
pip install -r requirements.txt
```

## Usage

### Simulate the process (Fe²⁺ and IP₆ over time)

```bash
python -m alaoline simulate --hours 36 --atmosphere n2
python -m alaoline simulate --hours 36 --atmosphere air
```

### Generate the 2⁽⁵⁻¹⁾ fractional-factorial design

```bash
python -m alaoline design --runs 16 --center 3 --out design.csv
```

### Fit a response surface from experimental data

```bash
python -m alaoline fit --data experiments.csv --out model.json
```

## As a library

```python
from alaoline.process import AlaoProcess
p = AlaoProcess(atmosphere="n2")
p.run(hours=36)
print(p.fe2_fraction, p.ip6_mg_100g)
```

## Testing

```bash
python -m pytest tests/ -q
```

## License

MIT. See [LICENSE](LICENSE).
