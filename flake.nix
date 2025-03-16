{
  description = "A Scalable Probabilistic Exact Model Counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    arjun = {
      url = "github:itepastra/arjun/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:itepastra/cryptominisat/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    sbva = {
      url = "github:itepastra/sbva/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    approxmc = {
      url = "github:itepastra/approxmc/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    breakid = {
      url = "github:itepastra/breakid/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      arjun,
      approxmc,
      breakid,
      cryptominisat,
      sbva,
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
          breakid,
          approxmc,
          python3,
          python3Packages,
        }:
        stdenv.mkDerivation {
          name = "ganak";
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
            breakid
            approxmc
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
            breakid = breakid.packages.${system}.breakid;
            approxmc = approxmc.packages.${system}.approxmc;
          };
        in
        {
          inherit ganak;
          default = ganak;
        }
      );
    };
}
