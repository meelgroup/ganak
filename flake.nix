{
  description = "A Scalable Probabilistic Exact Model Counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    # meel-cadiback is not used by ganak directly, but our deps
    # (cryptominisat, arjun, approxmc) pull it in. Declare it here as a single anchor
    # so the whole tree shares one instead of duplicating it in flake.lock.
    meel-cadiback = {
      url = "github:meelgroup/cadiback/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:msoos/cryptominisat/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.meel-cadiback.follows = "meel-cadiback";
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
      inputs.meel-cadiback.follows = "meel-cadiback";
      inputs.meel-cadical.follows = "meel-cadiback/meel-cadical";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
      inputs.treedecomp.follows = "treedecomp";
    };
    approxmc = {
      url = "github:meelgroup/approxmc/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.meel-cadiback.follows = "meel-cadiback";
      inputs.arjun.follows = "arjun";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
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
      meel-cadiback,
      arjun,
      approxmc,
      # breakid,
      cryptominisat,
      sbva,
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
