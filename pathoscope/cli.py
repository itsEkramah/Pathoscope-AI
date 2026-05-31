import sys
import argparse
from pathlib import Path
from loguru import logger
from pathoscope.pipeline import PipelineCoordinator
from pathoscope.utils.logger import setup_logger

def main():
    """
    Main entry point for the PathoScope AI CLI.
    Configures argument parsing, subcommands, and triggers execution.
    """
    parser = argparse.ArgumentParser(
        description="PathoScope AI: An Automated Viral Functional Genomics Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Global arguments
    parser.add_argument(
        "-c", "--config",
        type=str,
        default="config/default_config.yaml",
        help="Path to the YAML configuration file."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable extremely detailed debug-level logging."
    )
    
    subparsers = parser.add_subparsers(dest="command", required=True, help="Pipeline command to execute")
    
    # 'run' subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Execute the full end-to-end PathoScope AI pipeline."
    )
    run_parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to the raw viral sequence file (FASTA, FASTQ, CSV, TXT)."
    )
    run_parser.add_argument(
        "-o", "--outdir",
        type=str,
        default="results",
        help="Directory to write all versioned results directories."
    )
    
    # 'preprocess' subcommand
    pre_parser = subparsers.add_parser(
        "preprocess",
        help="Execute ONLY the sequence preprocessing and quality control step."
    )
    pre_parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to the raw sequence file."
    )
    pre_parser.add_argument(
        "-o", "--outdir",
        type=str,
        default="results",
        help="Directory to write preprocessed results."
    )

    args = parser.parse_args()

    # Configure global console logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(level=log_level)
    
    logger.info("Initializing PathoScope AI command line runner...")

    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Specified configuration file does not exist: {config_path}")
        sys.exit(1)

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Specified input sequence file does not exist: {input_path}")
        sys.exit(1)
        
    outdir_path = Path(args.outdir)

    try:
        coordinator = PipelineCoordinator(config_path=config_path, output_dir=outdir_path)
        
        if args.command == "run":
            logger.info("Triggering full pipeline execution workflow...")
            coordinator.execute_pipeline(input_path)
            
        elif args.command == "preprocess":
            logger.info("Triggering preprocessing-only workflow...")
            coordinator.initialize_workspace()
            try:
                coordinator.run_preprocessing(input_path)
                coordinator.write_run_metadata(input_path, status="SUCCESS")
            except Exception as e:
                coordinator.write_run_metadata(input_path, status="FAILED")
                raise e
            finally:
                coordinator.clean_temporary_files()
                
    except Exception as e:
        logger.exception(f"Execution failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
