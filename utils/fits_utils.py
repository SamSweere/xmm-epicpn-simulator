import numpy as np
from astropy.io import fits


def filter_evt_pattern_type(eventlist_path, max_pattern_type=4, verbose=True):
    # Load the eventlist
    with fits.open(eventlist_path) as hdu:
        # Load the event header and data
        events_header = hdu['EVENTS'].header
        events_data = np.array(hdu['EVENTS'].data)

        # import matplotlib.pyplot as plt
        # plt.hist(events_data['TYPE'])
        # plt.show()
        #
        # print("Type counts:", np.unique(events_data['TYPE'], return_counts=True))
        # print("Total len:", len(events_data['TYPE']))
        # print("Count 0", len(events_data['TYPE'] == 0))
        # print("Count 1:4", len(events_data['TYPE'] <= 4) - len(events_data['TYPE'] == 0))

        # Filter the events data, remove every event with pattern type > max pattern type
        filtered_events_data = events_data[events_data['TYPE'] <= max_pattern_type]

        # Since we filtered the events, set the patterns to 0
        for i in range(max_pattern_type + 1, 13):
            events_header.set(f'NGRAD{i}', 0)
            events_header.set(f'NPGRA{i}', 0)

        events_header.set('NAXIS2', len(filtered_events_data))
        # # Remove the checksum and datasum since these changed
        # events_header.remove('CHECKSUM')
        # events_header.remove('DATASUM')

        # Create new binary table from the filtered events
        filtered_events = fits.BinTableHDU(data=filtered_events_data, header=events_header)

        # Update the events with the new binary table
        hdu['EVENTS'] = filtered_events
        # Update history
        hdu['PRIMARY'].header['HISTORY'] = f"Removed all events with pattern type > {max_pattern_type}"

        # Overwrite the old eventlist
        hdu.writeto(eventlist_path, overwrite=True, checksum=True)

    return eventlist_path
