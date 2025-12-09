import sys
from casacore.tables import table

def remove_antennas(ms_path, antennas_to_remove):
    try:
        with table(f"{ms_path}/POINTING", readonly=False) as pointing_table:
            rows_to_remove = [
                i for i, antenna in enumerate(pointing_table)
                if antenna['NAME'] in antennas_to_remove
            ]
            if rows_to_remove:
                pointing_table.removerows(rows_to_remove)
                print(f"Removed rows from POINTING: {rows_to_remove}")

        with table(f"{ms_path}/ANTENNA", readonly=False) as antenna_table:
            rows_to_remove = [
                i for i, antenna in enumerate(antenna_table)
                if antenna['STATION'] in antennas_to_remove
            ]
            if rows_to_remove:
                antenna_table.removerows(rows_to_remove)
                print(f"Removed rows from ANTENNA: {rows_to_remove}")

        with table(f"{ms_path}/FLAG_CMD", readonly=False) as flag_table:
            rows_to_remove = [
                i for i, antenna in enumerate(flag_table)
                if antenna['STATION'] in antennas_to_remove
            ]
            if rows_to_remove:
                antenna_table.removerows(rows_to_remove)
                print(f"Removed rows from FLAG_CMD: {rows_to_remove}")

    except Exception as e:
        print(f"Error processing Measurement Set: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python antennakiller.py <MeasurementSetPath>")
        sys.exit(1)

    msin = sys.argv[1]
    antennas_to_remove = ['C07:31', 'S05:32', 'E01']
    remove_antennas(msin, antennas_to_remove)