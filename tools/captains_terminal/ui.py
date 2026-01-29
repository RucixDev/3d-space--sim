"""Tiny text UI helpers (stdlib only)."""

from __future__ import annotations

from typing import Iterable, List, Optional
import shutil
import textwrap


def term_width(default: int = 88) -> int:
    try:
        return max(60, min(120, shutil.get_terminal_size((default, 20)).columns))
    except Exception:
        return default


def rule(ch: str = "-", width: Optional[int] = None) -> str:
    w = term_width() if width is None else width
    return ch * max(0, w)


def header(title: str, *, width: Optional[int] = None) -> str:
    w = term_width() if width is None else width
    t = f" {title.strip()} "
    if len(t) >= w:
        return title.strip()
    pad = (w - len(t)) // 2
    return ("=" * pad) + t + ("=" * (w - len(t) - pad))


def wrap_lines(text: str, *, width: Optional[int] = None) -> List[str]:
    w = term_width() if width is None else width
    return textwrap.fill(text, width=w).splitlines()


def fmt_credits(n: int) -> str:
    return f"{int(n):,} cr"


def prompt_choice(prompt: str, *, allow_blank: bool = False) -> str:
    while True:
        try:
            s = input(prompt)
        except EOFError:
            return ""
        if s.strip() == "" and not allow_blank:
            continue
        return s.strip()


def prompt_int(prompt: str, *, min_v: int, max_v: int, allow_blank: bool = False) -> Optional[int]:
    while True:
        s = prompt_choice(prompt, allow_blank=allow_blank)
        if s == "" and allow_blank:
            return None
        try:
            v = int(s)
        except ValueError:
            print(f"Please enter an integer between {min_v} and {max_v}.")
            continue
        if v < min_v or v > max_v:
            print(f"Out of range. Enter {min_v}..{max_v}.")
            continue
        return v


def bullet_list(lines: Iterable[str], *, prefix: str = "- ") -> str:
    out: List[str] = []
    for ln in lines:
        out.append(prefix + ln)
    return "\n".join(out)
