// Headless test of the in-browser d-DNNF pipeline, run under Node.js.
//
// It exercises the exact same chain the page does, but without a browser:
//   ganak --compile (raw d4 .nnf)  ->  ddnnf-cleanup (strict d4)  ->  ddnnf2dot (DOT)
// using the very same ganak.js / ddnnf-cleanup.js / ddnnf2dot.js Emscripten
// modules that ship in this directory. Useful as a smoke test after rebuilding
// the wasm.
//
// Both .js files are non-MODULARIZE Emscripten CLI modules that expect a global
// `Module`. In the browser we give each its own Web Worker (own global scope);
// here we run each in a fresh function scope with `Module` already in scope,
// which is the equivalent trick (a plain `require` would let the module's own
// `var Module` shadow ours).
//
// Usage:
//   node chain.js [input.cnf] [out.dot]
//   node chain.js                      # uses the built-in pwmc example, prints DOT

const fs = require('fs');
const path = require('path');

const HERE = __dirname;

// Default example matches the page's textarea (a projected weighted CNF).
const DEFAULT_CNF = `c t pwmc
p cnf 3 2
c p weight 1 4/10 0
c p weight 2 6/10 0
c p weight 3 8/10 0
1 2 0
-3 -1 -2 0
c p show 1 2 3 0
`;

// Run one Emscripten CLI tool once: write inFiles into its in-memory FS, call
// main with args, return the contents of outFile (or null). Mirrors the worker
// in index.html.
function runTool(jsFile, args, inFiles, outFile) {
  return new Promise((resolve, reject) => {
    const code = fs.readFileSync(path.join(HERE, jsFile), 'utf8');
    const Module = {
      noInitialRun: true,
      arguments: [],
      print: t => process.stderr.write('[' + jsFile + '] ' + t + '\n'),
      printErr: t => {
        if (!t.includes('__syscall_getrusage')) process.stderr.write('[' + jsFile + '] ' + t + '\n');
      },
      onRuntimeInitialized: function () {
        try {
          for (const p in inFiles) Module.FS.writeFile(p, inFiles[p]);
          Module.callMain(args);
          resolve(outFile ? Module.FS.readFile(outFile, { encoding: 'utf8' }) : null);
        } catch (e) { reject(e); }
      }
    };
    // Keep the tool from auto-running main off our argv.
    process.argv = [process.argv[0]];
    new Function('Module', 'require', '__dirname', '__filename', code)(
      Module, require, HERE, path.join(HERE, jsFile));
  });
}

(async () => {
  const inPath = process.argv[2];
  const outPath = process.argv[3];
  const cnf = inPath ? fs.readFileSync(inPath, 'utf8') : DEFAULT_CNF;

  // 1. ganak --compile: dump the raw d4 .nnf (faithful, --weak 0 default).
  //    --mode 1 (weighted rational) handles both weighted and plain CNFs.
  const raw = await runTool('ganak.js',
    ['--verb', '1', '--mode', '1', '--compile', '/out.nnf', '/input.cnf'],
    { '/input.cnf': cnf }, '/out.nnf');
  process.stderr.write('\n=== RAW NNF (' + raw.length + ' bytes) ===\n');

  // 2. ddnnf-cleanup: drop dead nodes, renumber root=1 contiguous -> strict d4.
  const cleaned = await runTool('ddnnf-cleanup.js',
    ['/in.nnf', '/out.nnf'], { '/in.nnf': raw }, '/out.nnf');

  // 3. ddnnf2dot: render the cleaned circuit as Graphviz DOT.
  const dot = await runTool('ddnnf2dot.js',
    ['/in.nnf', '/out.dot'], { '/in.nnf': cleaned }, '/out.dot');

  if (outPath) {
    fs.writeFileSync(outPath, dot);
    process.stderr.write('\nwrote DOT to ' + outPath + '\n');
  } else {
    process.stdout.write(dot);
  }
  process.exit(0);
})().catch(e => { process.stderr.write('FAIL: ' + e + '\n'); process.exit(1); });
