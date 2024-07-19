import cdsapi

c = cdsapi.Client()

c.retrieve(
    'satellite-carbon-dioxide',
    {
        'processing_level': 'level_2',
        'sensor_and_algorithm': 'tanso2_fts2_srfp',
        'year': '2019',
        'month': '12',
        'day': '16',
        'version': '2.0.0',
        'format': 'zip',
    },
    'download.zip')
