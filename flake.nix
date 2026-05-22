{
  description = "A Scalable Probabilistic Exact Model Counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    # cadiback and evalmaxsat are not used by ganak directly, but our deps
    # (cryptominisat, arjun, approxmc) pull them in. Declare them here as
    # single anchors so the whole tree shares one of each instead of
    # duplicating them many times in flake.lock.
    cadiback = {
      url = "github:meelgroup/cadiback/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    evalmaxsat = {
      url = "github:meelgroup/EvalMaxSAT/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:msoos/cryptominisat/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.cadiback.follows = "cadiback";
    };
    sbva = {
      url = "github:meelgroup/sbva/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treedecomp = {
      url = "github:meelgroup/treedecomp/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    arjun = {
      url = "github:meelgroup/arjun/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.cadiback.follows = "cadiback";
      inputs.cadical.follows = "cadiback/cadical";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
      inputs.evalmaxsat.follows = "evalmaxsat";
      inputs.treedecomp.follows = "treedecomp";
    };
    approxmc = {
      url = "github:meelgroup/approxmc/master";
      inputs.nixpkgs.follows = "nixpkgs";
      # The cadiback follows is ignored on the old pre-dedup approxmc master
      # (which has no cadiback input) but activates once that flake declares
      # it. Everything else matches today.
      inputs.cadiback.follows = "cadiback";
      inputs.arjun.follows = "arjun";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
      inputs.evalmaxsat.follows = "evalmaxsat";
      inputs.treedecomp.follows = "treedecomp";
    };
    # breakid = {
    #   url = "github:meelgroup/breakid/master";
    #   inputs.nixpkgs.follows = "nixpkgs";
    # };
  };

  outputs =
    {
      self,
      nixpkgs,
      cadiback,
      arjun,
      approxmc,
      # breakid,
      cryptominisat,
      sbva,
      evalmaxsat,
      treedecomp,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});

      ganak-package =
        {
          stdenv,
          cmake,
          pkg-config,
          gmp,
          mpfr,
          flint3,
          zlib,
          autoPatchelfHook,
          cryptominisat,
          arjun,
          sbva,
          # breakid,
          approxmc,
          python3,
          python3Packages,
          # evalmaxsat,
          treedecomp,
        }:
        stdenv.mkDerivation {
          name = "ganak";
          cmakeFlags = [
            "-Dcryptominisat5_DIR=${cryptominisat}/lib/cmake/cryptominisat5"
            "-Darjun_DIR=${arjun}/lib/cmake/arjun"
            "-Dapproxmc_DIR=${approxmc}/lib/cmake/approxmc"
            "-Dtreedecomp_DIR=${treedecomp}/lib/cmake/treedecomp"
          ];
          nativeBuildInputs = [
            cmake
            pkg-config
            autoPatchelfHook
            python3
            python3Packages.numpy
          ];
          buildInputs = [
            gmp
            mpfr
            flint3
            zlib
            cryptominisat
            arjun
            sbva
            # breakid
            approxmc
            # evalmaxsat
            treedecomp
          ];
          src = ./.;
        };
    in
    {
      # checks

      devShells = forAllSystems (
        system:
        let
          pkgs = nixpkgsFor.${system};
        in
        {
          default = pkgs.mkShell {
            packages = [
              pkgs.cmake
              pkgs.pkg-config
              pkgs.gmp
            ];
          };
        }
      );

      packages = forAllSystems (
        system:
        let
          ganak = nixpkgsFor.${system}.callPackage ganak-package {
            cryptominisat = cryptominisat.packages.${system}.cryptominisat;
            arjun = arjun.packages.${system}.arjun;
            sbva = sbva.packages.${system}.sbva;
            # breakid = breakid.packages.${system}.breakid;
            approxmc = approxmc.packages.${system}.approxmc;
            # evalmaxsat = evalmaxsat.packages.${system}.evalmaxsat;
            treedecomp = treedecomp.packages.${system}.treedecomp;
          };
        in
        {
          inherit ganak;
          default = ganak;
        }
      );
    };
}
