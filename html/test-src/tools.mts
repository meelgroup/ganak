// Runs the Emscripten CLI tools in Node, one module per call, as the page's worker does.

import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import path from "node:path";
import { fileURLToPath } from "node:url";

export const HTML = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const req = createRequire(import.meta.url);

export interface Run {
  status: number;
  stdout: string;
  stderr: string;
  files: Record<string, Uint8Array>;
}

// One Emscripten CLI module per call, as the page's worker does.
export function runTool(tool: string, args: string[], files: Record<string, string | Uint8Array>,
    outFiles: string[] = []): Promise<Run> {
  return new Promise((resolve) => {
    const js = path.join(HTML, tool + ".js");
    const out: string[] = [];
    const err: string[] = [];
    const mod: any = {
      noInitialRun: true,
      arguments: [],
      print: (t: string) => out.push(t),
      printErr: (t: string) => err.push(t),
      quit: (_s: number, e: unknown) => { throw e; },
      onRuntimeInitialized: () => {
        for (const [p, d] of Object.entries(files)) mod.FS.writeFile(p, d);
        let status = 0;
        try {
          const r = mod.callMain([...args]);  // callMain changes the array
          if (typeof r === "number") status = r;
        } catch (x: any) {
          if (x && typeof x.status === "number") status = x.status;
          else { err.push(String(x)); status = -1; }
        }
        process.exitCode = 0;  // Emscripten sets it on exit(); the tests decide
        const got: Record<string, Uint8Array> = {};
        for (const p of outFiles) { try { got[p] = mod.FS.readFile(p); } catch { /* not written */ } }
        resolve({ status, stdout: out.join("\n"), stderr: err.join("\n"), files: got });
      },
    };
    const argv = process.argv;
    process.argv = [argv[0]];
    new Function("Module", "require", "__dirname", "__filename", readFileSync(js, "utf8"))(mod, req, HTML, js);
    process.argv = argv;
  });
}

