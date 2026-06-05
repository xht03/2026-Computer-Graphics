import os
import subprocess
import tempfile

BLENDER_PATH = os.getenv("BLENDER_PATH", "blender")


def _parse_result(returncode: int, stdout: str, stderr: str) -> dict:
    """Determine true success: Blender exits 0 even on Python script errors,
    so we check for our sentinel string and absence of Traceback."""
    has_sentinel = "ALL_RENDERS_COMPLETE" in stdout
    has_traceback = "Traceback (most recent call last)" in stdout
    success = returncode == 0 and has_sentinel and not has_traceback

    error_detail = ""
    if has_traceback:
        lines = stdout.splitlines()
        tb_start = next((i for i, l in enumerate(lines) if "Traceback" in l), None)
        if tb_start is not None:
            error_detail = "\n".join(lines[tb_start:tb_start + 20])

    return {
        "success": success,
        "stdout": stdout,
        "stderr": stderr,
        "returncode": returncode,
        "error_detail": error_detail,
    }


def run_script(script_code: str, timeout: int = 180) -> dict:
    """Execute a bpy script string in Blender headless mode."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False, encoding="utf-8") as f:
        f.write(script_code)
        script_path = f.name

    try:
        cmd = [BLENDER_PATH, "--background", "--python", script_path]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, env={**os.environ}
        )
        return _parse_result(result.returncode, result.stdout, result.stderr)
    except subprocess.TimeoutExpired:
        return {"success": False, "stdout": "", "stderr": f"Timed out after {timeout}s",
                "returncode": -1, "error_detail": ""}
    except FileNotFoundError:
        return {"success": False, "stdout": "", "stderr": f"Blender not found at '{BLENDER_PATH}'.",
                "returncode": -1, "error_detail": ""}
    finally:
        os.unlink(script_path)


def run_script_file(script_path: str, timeout: int = 180) -> dict:
    """Execute an existing bpy script file in Blender headless mode."""
    try:
        cmd = [BLENDER_PATH, "--background", "--python", script_path]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, env={**os.environ}
        )
        return _parse_result(result.returncode, result.stdout, result.stderr)
    except subprocess.TimeoutExpired:
        return {"success": False, "stdout": "", "stderr": f"Timed out after {timeout}s",
                "returncode": -1, "error_detail": ""}
    except FileNotFoundError:
        return {"success": False, "stdout": "", "stderr": f"Blender not found at '{BLENDER_PATH}'.",
                "returncode": -1, "error_detail": ""}
