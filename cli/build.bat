@echo off
REM Windows build. Run from a "x64 Native Tools Command Prompt for VS" (MSVC)
REM or from an MSYS2 / MinGW-w64 shell (gcc). No zlib needed: this build
REM reads plain text only, which is what --format plink2 output is.
where cl >nul 2>nul
if %errorlevel%==0 (
  cl /O2 /DPHYSMERGE_NO_ZLIB /Fephysmerge.exe physmerge.c
) else (
  gcc -O2 -std=c99 -DPHYSMERGE_NO_ZLIB -o physmerge.exe physmerge.c
)
