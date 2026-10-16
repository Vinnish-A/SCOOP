import pytest
from scoop.atlas import scvi_training_options


def test_short_schedule_native_options_and_cli_precedence():
    config = {'atlas': {'scvi': {'max_epochs': 60, 'batch_size': 2048,
                                'train_kwargs': {'datasplitter_kwargs': {'num_workers': 2}, 'accelerator': 'cpu'}}}}
    options = scvi_training_options(config)
    assert options['plan_kwargs']['n_epochs_kl_warmup'] == 20
    assert options['batch_size'] == 2048
    assert options['datasplitter_kwargs']['num_workers'] == 2
    assert options['accelerator'] == 'cpu'
    assert scvi_training_options(config, batch_size=512)['batch_size'] == 512
    config['atlas']['scvi']['plan_kwargs'] = {'n_epochs_kl_warmup': 400}
    with pytest.raises(ValueError, match='KL warmup'):
        scvi_training_options(config)
    options = scvi_training_options(config, kl_warmup_epochs=10)
    assert options['plan_kwargs']['n_epochs_kl_warmup'] == 10
