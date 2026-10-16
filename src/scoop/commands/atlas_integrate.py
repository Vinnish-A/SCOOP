from __future__ import annotations

import argparse

from scoop.atlas import integrate


def main() -> int:
    parser = argparse.ArgumentParser(description="Run scVI integration and clustering in X_scVI.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--max-epochs", type=int)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--accelerator")
    parser.add_argument("--devices", help="Number of devices or auto; use train_kwargs.devices for a device list")
    parser.add_argument("--kl-warmup-epochs", type=int)
    args = parser.parse_args()
    print(
        integrate(config_path=args.config, input_path=args.input, output_path=args.output, max_epochs=args.max_epochs,
                  batch_size=args.batch_size, accelerator=args.accelerator,
                  devices=int(args.devices) if args.devices and args.devices.isdigit() else args.devices,
                  kl_warmup_epochs=args.kl_warmup_epochs)
    )
    return 0
