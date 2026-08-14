# Getting Started

Welcome! This repository is your workspace for the course. All programming happens inside GitHub Codespaces, a cloud development environment that comes with every tool we need already installed. You do not have to install anything on your own computer; a web browser and a GitHub account are enough.

## Opening the Repository in a Codespace

Sign in to GitHub and navigate to this repository. Click the green button labeled "Code" near the top of the page, select the "Codespaces" tab, and click "Create codespace on main". GitHub then builds your personal copy of the course environment. The first launch takes a few minutes while the container image is downloaded, so this is a good moment to get a coffee. Once it finishes, VS Code opens directly in your browser and you are ready to work.

A codespace is persistent. The next time you return to the repository, your existing codespace appears under the same "Codespaces" tab, and you should reopen it instead of creating a new one. Creating many codespaces wastes your free monthly quota, and your files live inside the codespace, so reopening the same one means all your work is exactly where you left it.

## Using Codespaces from VS Code on Your Computer

If you prefer working in a locally installed editor over the browser, you can connect VS Code on your own machine to your codespace. Install VS Code from https://code.visualstudio.com, open the Extensions view, and install the extension called "GitHub Codespaces". Sign in to GitHub when the extension asks you to. Then open the Command Palette with Cmd+Shift+P on macOS or Ctrl+Shift+P on Windows and Linux, run the command "Codespaces: Connect to Codespace", and pick the codespace you created earlier. The editor now runs on your computer while your code still compiles and runs in the cloud container, so everything else works exactly as in the browser.

## What Is Inside the Environment

The container is based on Ubuntu and ships the LLVM 22 toolchain, so `clang`, `clang++`, `clangd`, `clang-format` and the `lldb` debugger are all available under their usual names, together with CMake. The SQLite library and its command line tool `sqlite3` are also installed; the introduction uses them as an example of a C library. We will compile everything with `-std=c++23`.

## Checking Your Setup

To confirm that your environment is set up correctly, open a terminal inside VS Code (menu Terminal, then New Terminal) and run the command `getpasscode`. It prints a passcode, you setup everything perfectly. If the terminal reports that the command is not found, you are not inside the course container — make sure you opened the repository in a codespace as described above.

## Where to Start

Once VS Code is open, in the browser or locally, start by working through [introduction-to-programming-in-cpp.md](https://github.com/gtjusila/cio-ss-2026/blob/main/introduction-to-programming-in-cpp.md). It assumes you have programmed in Python before and introduces the parts of C++ we use in this course.
