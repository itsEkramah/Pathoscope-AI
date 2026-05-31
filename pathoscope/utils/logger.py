import sys
from pathlib import Path
from loguru import logger

def setup_logger(log_file: Path = None, level: str = "INFO"):
    """
    Configure loguru logger for console and optional file outputs.
    Removes default handlers and defines custom structured output layout.
    """
    # Remove standard default logger output to prevent duplicate formats
    logger.remove()
    
    # Custom colored format for console printout
    console_format = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
        "<level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
        "<level>{message}</level>"
    )
    
    logger.add(sys.stderr, format=console_format, level=level, colorize=True)
    
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        # Custom file format (including datetime and process info)
        file_format = (
            "{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | "
            "{name}:{function}:{line} - {message}"
        )
        logger.add(
            log_file,
            format=file_format,
            level=level,
            rotation="10 MB",
            compression="zip",
            encoding="utf-8"
        )
        logger.info(f"Log file initialized at: {log_file}")
    
    return logger
