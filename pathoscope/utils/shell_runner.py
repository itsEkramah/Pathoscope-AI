import subprocess
import sys
from pathlib import Path
from typing import List, Union
from loguru import logger

class ShellRunnerError(Exception):
    """Base exception for shell runner failures."""
    pass

class SubprocessExecutionError(ShellRunnerError):
    """Raised when an external subprocess execution returns a non-zero exit code."""
    def __init__(self, command: str, exit_code: int, stdout: str, stderr: str):
        self.command = command
        self.exit_code = exit_code
        self.stdout = stdout
        self.stderr = stderr
        super().__init__(
            f"Command failed with exit code {exit_code}.\n"
            f"Command: {command}\n"
            f"Stderr: {stderr.strip()[:1000]}"
        )


def execute_cmd(
    cmd: Union[str, List[str]],
    cwd: Path = None,
    env: dict = None
) -> subprocess.CompletedProcess:
    """
    Safely executes a shell command using subprocess.run().
    Logs the command, captures outputs, and raises SubprocessExecutionError on failure.
    """
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    logger.debug(f"Executing subprocess command: {cmd_str}")
    
    try:
        # Run command with shell=True if passed as string, otherwise shell=False
        is_shell = isinstance(cmd, str)
        
        result = subprocess.run(
            cmd,
            shell=is_shell,
            cwd=str(cwd) if cwd else None,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace"
        )
        
        if result.returncode != 0:
            logger.error(f"Command failed (exit code: {result.returncode})")
            logger.error(f"Stderr: {result.stderr.strip()}")
            raise SubprocessExecutionError(
                command=cmd_str,
                exit_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr
            )
            
        logger.debug(f"Command completed successfully (exit code 0)")
        return result
        
    except subprocess.SubprocessError as se:
        logger.error(f"Subprocess system failure for command '{cmd_str}': {se}")
        raise ShellRunnerError(f"System failed to execute command: {se}")
