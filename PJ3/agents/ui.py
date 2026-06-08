"""Centralized terminal UI (rich-based).

All user-facing output goes through this module so the CLI has one consistent,
Claude-Code-like look: a shared Console, status spinners for long steps, colored
status lines, a rounded review panel, and a clean action menu. Agents print their
incidental status via ui.log() so it cooperates with the live spinners (builtin
print() would corrupt the live region).
"""
from __future__ import annotations

import os
from contextlib import contextmanager

from rich.console import Console, Group
from rich.panel import Panel
from rich.rule import Rule
from rich.text import Text

console = Console(highlight=False)

_VIEWS = ("front", "side", "top", "iso")
_SEV_COLOR = {"error": "red", "warning": "yellow", "info": "cyan"}
_DIFF_COLOR = {"ADD": "green", "REMOVE": "red", "CHANGE": "yellow"}


# ── low-level helpers ─────────────────────────────────────────────────────────

def log(text: str) -> None:
    """Incidental/secondary status from agents (dim). Safe under a live spinner."""
    console.print(f"  [dim]{text}[/]")


def ok(label: str, detail: str = "") -> None:
    line = f"[green]✓[/] {label}"
    if detail:
        line += f" [dim]· {detail}[/]"
    console.print(line)


def warn(label: str, detail: str = "") -> None:
    line = f"[yellow]●[/] {label}"
    if detail:
        line += f" [dim]· {detail}[/]"
    console.print(line)


def err(label: str, detail: str = "") -> None:
    line = f"[red]✗[/] {label}"
    if detail:
        line += f" [dim]· {detail}[/]"
    console.print(line)


def note(text: str) -> None:
    console.print(f"  [dim]{text}[/]")


@contextmanager
def working(label: str, spinner: bool = True):
    """Wrap a long, blocking step. Shows a spinner unless disabled (e.g. --verbose,
    where the agent's own debug output would fight the live region)."""
    if spinner:
        with console.status(f"[bold cyan]{label}[/]…", spinner="dots"):
            yield
    else:
        console.print(f"  [dim]▸ {label}…[/]")
        yield


# ── composite views ───────────────────────────────────────────────────────────

def banner(description: str, *, vlm: bool, geom: bool, interactive: bool,
           max_iterations: int) -> None:
    def flag(name, on):
        return f"[green]{name} ✓[/]" if on else f"[dim]{name} ✗[/]"

    mode = "interactive" if interactive else f"auto · max {max_iterations}"
    console.print()
    console.print(f"  [bold]◆ PJ3 Blender Agent[/]")
    console.print(f"  [bold cyan]{description}[/]")
    console.print(f"  {flag('vlm', vlm)}   {flag('geom', geom)}   [dim]· {mode}[/]")
    console.print(Rule(style="grey30"))


def _issue_text(i: dict) -> Text:
    sev = i.get("severity", "info")
    color = _SEV_COLOR.get(sev, "white")
    kind = i.get("kind", "parameter")
    src = i.get("source", "?")
    t = Text("  ")
    t.append("● ", style=color)
    t.append(f"{sev:<7}", style=color)
    t.append(f" {src}/{kind} ", style="dim")
    t.append(" " + i.get("message", ""), style="default")
    return t


def review(iteration: int, iter_dir: str, issues: list[dict]) -> None:
    """Render the per-round review panel: renders + detected issues."""
    body: list = []

    have = [v for v in _VIEWS if os.path.exists(os.path.join(iter_dir, f"{v}.png"))]
    rend = Text("  ")
    rend.append("Renders  ", style="bold")
    if have:
        rend.append("  ".join(have), style="default")
        rend.append(f"   ·  {iter_dir}/", style="dim")
    else:
        rend.append("none — the script crashed", style="red")
    body.append(rend)
    body.append(Text(""))

    head = Text("  ")
    head.append("Issues  ", style="bold")
    if issues:
        head.append(str(len(issues)), style="default")
    else:
        head.append("none", style="green")
    body.append(head)

    order = {"error": 0, "warning": 1, "info": 2}
    for i in sorted(issues, key=lambda x: order.get(x.get("severity", "info"), 2)):
        body.append(_issue_text(i))

    title = f"[bold]Round {iteration}[/] [dim]· {os.path.basename(iter_dir)}[/]"
    console.print()
    console.print(Panel(Group(*body), title=title, title_align="left",
                        border_style="grey37", padding=(0, 1)))


def menu() -> None:
    console.print(
        "  [bold cyan]a[/] 接受并完成    "
        "[bold cyan]f[/] 修复检测到的问题    "
        "[bold cyan]c[/] 提出修改意见    "
        "[bold cyan]q[/] 退出"
    )


def ask(prompt: str = "› ") -> str:
    try:
        return console.input(f"  [bold cyan]{prompt}[/]").strip()
    except EOFError:
        return ""


def plan_diff(changes: str) -> None:
    console.print("  [bold]Plan diff[/] [dim]→ applying incrementally[/]")
    for ln in changes.splitlines():
        head = ln.strip().split(" ", 1)[0]
        color = _DIFF_COLOR.get(head)
        console.print(f"    {ln}", style=color or "dim")


def footer(run_dir: str, success: bool) -> None:
    console.print(Rule(style="grey30"))
    if success:
        console.print(f"  [green]✓ Done[/] [dim]· {run_dir}[/]")
    else:
        console.print(f"  [red]✗ Stopped[/] [dim]· {run_dir}[/]")
    console.print()
