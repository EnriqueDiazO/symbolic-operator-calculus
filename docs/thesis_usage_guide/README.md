# Guía visual de uso para la tesis

## Objetivo

Esta guía acompaña el estudio simbólico del problema de Haseman con shift para
el caso actualmente implementado, $m=2$. Presenta el AST no conmutativo, los
kernels Fourier normalizados, las factorizaciones Wiener–Hopf relativas en las
dos direcciones, sus tres evidencias, la primera reducción de Schur, los kernels
laterales y combinado, y las proyecciones LaTeX.

El notebook fue generado y ejecutado sobre el commit
`3216465f6afd0fa97c627c73ba4fba95110dc179`.

## Archivos principales

- Notebook ejecutado:
  `notebooks/thesis_symbolic_operator_calculus_guide.ipynb`.
- Exportación HTML:
  `docs/thesis_usage_guide/thesis_symbolic_operator_calculus_guide.html`.
- Capturas reales: `docs/thesis_usage_guide/screenshots/`.

## Abrir y ejecutar

Desde la raíz del repositorio, puede abrirse interactivamente con:

```bash
PYENV_VERSION=3.10.14 PYTHONPATH=src jupyter lab
```

Para ejecutar todas las celdas desde un kernel limpio:

```bash
PYTHONDONTWRITEBYTECODE=1 \
PYENV_VERSION=3.10.14 \
PYTHONPATH="$PWD/src" \
python -m jupyter nbconvert \
  --to notebook \
  --execute \
  --inplace \
  --ExecutePreprocessor.timeout=180 \
  notebooks/thesis_symbolic_operator_calculus_guide.ipynb
```

El comando no usa `allow_errors`: cualquier celda fallida detiene la ejecución.

## Regenerar el HTML

Después de ejecutar el notebook y de disponer de las capturas:

```bash
PYENV_VERSION=3.10.14 \
PYTHONPATH="$PWD/src" \
python -m jupyter nbconvert \
  --to html \
  --HTMLExporter.embed_images=True \
  --output-dir docs/thesis_usage_guide \
  notebooks/thesis_symbolic_operator_calculus_guide.ipynb
```

La opción de exportación incorpora la galería al HTML autónomo. El notebook
mantiene referencias relativas a los PNG y no contiene imágenes base64.

Para una regeneración completa desde cero, el orden es:

1. ejecutar el notebook;
2. exportar una primera vez a HTML sin `--HTMLExporter.embed_images=True`;
3. generar las diez capturas con el comando de la sección siguiente;
4. repetir la exportación anterior con `--HTMLExporter.embed_images=True` para
   obtener el HTML final con la galería incorporada.

## Regenerar las capturas

Las capturas se obtuvieron del HTML ejecutado con Pyppeteer 2.0.0 y el Chrome
ya disponible en `/usr/bin/google-chrome`. El siguiente comando espera a que
MathJax termine, desplaza cada ancla real y captura el viewport correspondiente:

```bash
PYENV_VERSION=3.10.14 python - <<'PY'
import asyncio
from pathlib import Path

from pyppeteer import launch

captures = (
    ("01_portada_y_alcance.png", "sec-portada"),
    ("02_ast_no_conmutativo.png", "sec-ast"),
    ("03_kernels_fourier.png", "sec-fourier"),
    ("04_factorizacion_relativa_1_2.png", "sec-factor12"),
    ("05_acciones_relativas.png", "sec-actions-evidence"),
    ("06_correspondencias_fourier.png", "sec-fourier-evidence"),
    ("07_factorizacion_relativa_2_1.png", "sec-factor21"),
    ("08_reduccion_schur_m2.png", "sec-schur"),
    ("09_kernel_combinado_c22.png", "sec-c22"),
    ("10_renderizado_latex.png", "sec-latex"),
)


async def main():
    browser = await launch(
        executablePath="/usr/bin/google-chrome",
        headless=True,
        args=[
            "--no-sandbox",
            "--disable-gpu",
            "--hide-scrollbars",
            "--allow-file-access-from-files",
        ],
    )
    page = await browser.newPage()
    await page.setViewport(
        {"width": 1400, "height": 1050, "deviceScaleFactor": 1}
    )
    html = Path(
        "docs/thesis_usage_guide/"
        "thesis_symbolic_operator_calculus_guide.html"
    ).resolve().as_uri()
    await page.goto(html, {"waitUntil": "networkidle0", "timeout": 30000})
    await page.evaluate(
        """async () => {
            if (window.MathJax && MathJax.startup && MathJax.startup.promise) {
                await MathJax.startup.promise;
            }
        }"""
    )
    output = Path("docs/thesis_usage_guide/screenshots")
    for filename, anchor in captures:
        found = await page.evaluate(
            f"Boolean(document.getElementById('{anchor}'))"
        )
        if not found:
            raise RuntimeError(f"Ancla ausente: {anchor}")
        await page.evaluate(
            f"document.getElementById('{anchor}').scrollIntoView("
            "{block: 'start'})"
        )
        await asyncio.sleep(0.35)
        await page.screenshot(
            {"path": str(output / filename), "fullPage": False}
        )
    await browser.close()


asyncio.run(main())
PY
```

No se instalaron herramientas para esta fase y no se usaron mockups ni HTML
fabricado.

## Inventario de capturas

| PNG | Sección del notebook |
| --- | --- |
| `01_portada_y_alcance.png` | Portada y alcance actual |
| `02_ast_no_conmutativo.png` | Primer contacto con el AST |
| `03_kernels_fourier.png` | Kernels Fourier normalizados |
| `04_factorizacion_relativa_1_2.png` | Factorización relativa $1\to2$ |
| `05_acciones_relativas.png` | Evidencia de acciones |
| `06_correspondencias_fourier.png` | Evidencia Fourier independiente |
| `07_factorizacion_relativa_2_1.png` | Dirección $2\to1$ |
| `08_reduccion_schur_m2.png` | Primera reducción de Schur |
| `09_kernel_combinado_c22.png` | Kernels laterales y combinado |
| `10_renderizado_latex.png` | Renderizado y exportación LaTeX |

## Limitaciones

- El dominio Wiener–Hopf relativo y el AST de Schur siguen separados.
- $R_{1,1}$ es un regularizador formal, no un inverso computado.
- La equivalencia módulo compactos es metadato matemático validado, no una
  decisión automática de compacidad.
- No hay reducción general para $m>2$, inversión real de operadores ni
  evaluación numérica general.
- La guía demuestra el comportamiento del prototipo; no reemplaza las pruebas
  analíticas de la tesis.
