# Pathoscope-AI

Pathoscope-AI is an automated viral functional genomics and gene expression analysis pipeline designed for research, teaching, and clinical bioinformatics exploration. It supports raw sequence preprocessing, ORF prediction, functional annotation, pathway mapping, enrichment analysis, scientific visualization, and AI-assisted interpretation.

## Key Features

- Preprocesses FASTA/FASTQ input and quality filters viral sequences
- Predicts open reading frames (ORFs) and translates them into protein sequences
- Performs similarity-based functional annotation using local reference databases
- Maps proteins to pathways and Pfam domains
- Runs enrichment analysis and generates publication-ready visualizations
- Supports AI interpretation via OpenAI/Gemini with offline fallback
- Provides both CLI and Streamlit web interface

## Project Structure

- `pathoscope/` — main application package
- `config/default_config.yaml` — default pipeline configuration

- `requirements.txt` — Python dependency list
- `README.md` — this file

## Setup

1. Open PowerShell in the project root:

```powershell
cd "d:\PathoScope AI\Pathoscope-AI"
```

2. Activate the local virtual environment:

```powershell
.\.venv\Scripts\Activate.ps1
```

3. Install dependencies:

```powershell
python -m pip install -r requirements.txt
```

## Optional Configuration

The project can read `.env` from the repository root. Create or update `.env` with your API keys, for example:

```env
OPENAI_API_KEY=sk-...your-openai-key...
GEMINI_API_KEY=AIzaSy...your-gemini-key...
```

> Do not commit `.env` to version control.

## CLI Usage

### Run the full pipeline

```powershell
python -m pathoscope.cli run -i sample.fasta -o tmp_outdir
```

This creates a versioned run folder under `tmp_outdir`, including:

- `preprocessed/`
- `orfs/`
- `annotations/`
- `pathways/`
- `enrichment/`
- `visualizations/`
- `final_report/`
- `metadata.json`

### Run only preprocessing

```powershell
python -m pathoscope.cli preprocess -i sample.fasta -o tmp_outdir
```

## Streamlit App

Launch the web interface with:

```powershell
python -m streamlit run pathoscope/app.py
```

Then open the browser at:

- `http://localhost:8501`

## Notes

- The AI interpretation stage uses `OPENAI_API_KEY` by default.
- If no API key is available, the app falls back to offline rule-based interpretation.
- The pathway mapping stage relies on `data/reference/Pfam-A.hmm` and `data/reference/pathways_cache.db`.
- PDF export requires additional WeasyPrint system dependencies; HTML report generation works without it.

## Troubleshooting

- If a dependency import fails, verify you are using the `.venv` Python interpreter.
- If `loguru` is missing, re-run:

```powershell
python -m pip install loguru
```

- If Streamlit does not start on `8501`, try a different port:

```powershell
python -m streamlit run pathoscope/app.py --server.port 8502
```

## License

Pathoscope-AI is distributed under the MIT License.

