import os
import pickle
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np # Add numpy import
from datetime import datetime, timedelta
from pyubx2 import UBXReader
from collections import defaultdict
import math

# Define the ESF-RAW sensor data type map based on uBlox documentation (e.g., Table 15)
# Format: data_type_id: {"name": "field_name", "unit": "unit_str", "scale": scale_factor}
# All data values are 24-bit signed integers.
ESF_RAW_SENSOR_DEFINITION_MAP = {
    # Gyroscope (deg/s, scale 2^-12)
    0x05: {"name": "raw_gyro_z", "unit": "deg/s", "scale": 2**-12}, # Z-axis gyroscope
    0x0E: {"name": "raw_gyro_x", "unit": "deg/s", "scale": 2**-12}, # X-axis gyroscope
    0x0D: {"name": "raw_gyro_y", "unit": "deg/s", "scale": 2**-12}, # Y-axis gyroscope
    # Accelerometer (m/s^2, scale 2^-10)
    0x10: {"name": "raw_accel_x", "unit": "m/s^2", "scale": 2**-10}, # X-axis accelerometer
    0x11: {"name": "raw_accel_y", "unit": "m/s^2", "scale": 2**-10}, # Y-axis accelerometer
    0x12: {"name": "raw_accel_z", "unit": "m/s^2", "scale": 2**-10}, # Z-axis accelerometer
    # Temperature (deg C, scale 1e-2)
    0x0C: {"name": "raw_temp", "unit": "deg C", "scale": 1e-2},      # Temperature
    # Magnetometer (LSB, scale 1.0 - as per Table 15 not providing specific scaling for ESF-RAW)
    0x00: {"name": "raw_mag_x", "unit": "LSB", "scale": 1.0},       # X-axis magnetometer
    0x01: {"name": "raw_mag_y", "unit": "LSB", "scale": 1.0},       # Y-axis magnetometer
    0x02: {"name": "raw_mag_z", "unit": "LSB", "scale": 1.0},       # Z-axis magnetometer
}

def parse_ubx_file(filename):
    """
    Parses a UBX file and extracts position, satellite, and ESF data.
    """
    positions = []
    satellite_data_list = [] # Renamed to avoid conflict with pandas DataFrame
    esf_data_list = []       # Renamed to avoid conflict with pandas DataFrame
    esf_cal_list = []
    esf_alg_list = []
    esf_raw_list = []
    esf_ins_list = [] # Add list for ESF-INS data
    nav_att_list = [] # Add list for NAV-ATT data
    nav_dop_list = [] # Add list for NAV-DOP data
    esf_raw_stag_list = []   # List to temporarily store raw ESF-RAW data with sTag values
    raw_msgs = defaultdict(list)
    last_known_pvt_time = None
    first_pvt_time = None
    last_pvt_time = None

    print(f"Parsing UBX file: {filename}")
    try:
        with open(filename, "rb") as f:
            ubr = UBXReader(f)
            for raw_data, parsed_msg in ubr:
                msg_type = parsed_msg.identity
                raw_msgs[msg_type].append(parsed_msg)

                if msg_type == "NAV-PVT":
                    # print(parsed_msg)
                    if all(hasattr(parsed_msg, attr) for attr in ['lat', 'lon', 'year', 'month', 'day', 'hour', 'min', 'second']):
                        # Create base datetime object from YMDHMS
                        base_dt = datetime(
                            parsed_msg.year,
                            parsed_msg.month,
                            parsed_msg.day,
                            parsed_msg.hour,
                            parsed_msg.min,
                            parsed_msg.second
                        )
                        # Get nanoseconds, default to 0 if not present
                        nanos = getattr(parsed_msg, 'nano', 0)
                        # Add nanoseconds as a timedelta (timedelta takes microseconds)
                        # This correctly handles positive and negative nanosecond offsets
                        timestamp = base_dt + timedelta(microseconds=nanos // 1000)
                        
                        last_known_pvt_time = timestamp
                        # Track first and last PVT times for interpolation
                        if first_pvt_time is None:
                            first_pvt_time = timestamp
                        last_pvt_time = timestamp
                        
                        positions.append({
                            "time": timestamp,
                            "gnssFixOk": parsed_msg.gnssFixOk,
                            "heading": parsed_msg.headVeh, # This is from NAV-PVT, distinct from NAV-ATT heading
                            "speed": parsed_msg.gSpeed / 1000.0 * 3.6, # mm/s to km/h
                            "velN": parsed_msg.velN,
                            "velE": parsed_msg.velE,
                            "velD": parsed_msg.velD,
                            "lat": parsed_msg.lat,
                            "lon": parsed_msg.lon,
                            "hMSL": getattr(parsed_msg, 'hMSL', np.nan) / 1000.0,  # mm to m
                            "hAcc": getattr(parsed_msg, 'hAcc', np.nan) / 1000.0,  # mm to m
                            "vAcc": getattr(parsed_msg, 'vAcc', np.nan) / 1000.0,   # mm to m
                            "fixType": getattr(parsed_msg, 'fixType', -1) # Add fixType, default to -1 if not present
                        })

                elif msg_type == "NAV-DOP":
                    # Handle NAV-DOP message
                    current_time_for_msg = last_known_pvt_time
                    print(f'current_time_for_msg: {current_time_for_msg}')
                    if current_time_for_msg and all(hasattr(parsed_msg, attr) for attr in ['iTOW', 'gDOP', 'pDOP', 'tDOP', 'vDOP', 'hDOP', 'nDOP', 'eDOP']):
                        print(f"[DEBUG] NAV-DOP message at {current_time_for_msg} with iTOW: {parsed_msg.iTOW}")
                        # DOP values are scaled by 100 in UBX protocol, so we divide by 100 to get actual values
                        nav_dop_list.append({
                            "time": current_time_for_msg,
                            "iTOW": parsed_msg.iTOW,
                            "gDOP": parsed_msg.gDOP / 100.0,  # Geometric DOP
                            "pDOP": parsed_msg.pDOP / 100.0,  # Position DOP
                            "tDOP": parsed_msg.tDOP / 100.0,  # Time DOP
                            "vDOP": parsed_msg.vDOP / 100.0,  # Vertical DOP
                            "hDOP": parsed_msg.hDOP / 100.0,  # Horizontal DOP
                            "nDOP": parsed_msg.nDOP / 100.0,  # Northing DOP
                            "eDOP": parsed_msg.eDOP / 100.0   # Easting DOP
                        })

                elif msg_type == "NAV-ATT":
                    current_time_for_msg = last_known_pvt_time # Approximate time
                    # pyubx2 already applies the 1e-5 scaling for roll, pitch, heading and accuracies
                    if current_time_for_msg and all(hasattr(parsed_msg, attr) for attr in ['iTOW', 'roll', 'pitch', 'heading']):
                        nav_att_list.append({
                            "time": current_time_for_msg,
                            "iTOW": parsed_msg.iTOW,
                            "roll": parsed_msg.roll,
                            "pitch": parsed_msg.pitch,
                            "heading": parsed_msg.heading, # This is the heading from NAV-ATT
                            "accRoll": getattr(parsed_msg, 'accRoll', np.nan),
                            "accPitch": getattr(parsed_msg, 'accPitch', np.nan),
                            "accHeading": getattr(parsed_msg, 'accHeading', np.nan)
                        })
                    # else:
                        # print(f"[DEBUG] Skipping NAV-ATT (iTOW: {getattr(parsed_msg, 'iTOW', 'N/A')}) due to missing time or attributes.")

                elif msg_type == "NAV-SAT":
                    current_time_for_msg = last_known_pvt_time
                    if current_time_for_msg and hasattr(parsed_msg, 'numSvs'):
                        num_svs = parsed_msg.numSvs
                        sats_details = []
                        for i in range(1, num_svs + 1):
                            svid_attr = f"svId_{i:02d}"
                            cno_attr = f"cno_{i:02d}"
                            elev_attr = f"elev_{i:02d}"
                            azim_attr = f"azim_{i:02d}"
                            
                            if hasattr(parsed_msg, svid_attr) and hasattr(parsed_msg, cno_attr):
                                sats_details.append({
                                    "svId": getattr(parsed_msg, svid_attr),
                                    "cno": getattr(parsed_msg, cno_attr),
                                    "elev": getattr(parsed_msg, elev_attr, np.nan),
                                    "azim": getattr(parsed_msg, azim_attr, np.nan),
                                })
                        
                        satellite_data_list.append({
                            "time": current_time_for_msg,
                            "iTOW": parsed_msg.iTOW,
                            "numSvs": num_svs,
                            "sats_details": sats_details
                        })

                elif msg_type == "ESF-MEAS":
                    current_time_for_msg = last_known_pvt_time
                    if current_time_for_msg and hasattr(parsed_msg, 'numMeas') and parsed_msg.numMeas > 0:
                        measurements = {}
                        has_imu_data = False
                        for i in range(1, parsed_msg.numMeas + 1):
                            df_attr = f"dataField_{i:02d}"
                            dt_attr = f"dataType_{i:02d}"
                            
                            if hasattr(parsed_msg, df_attr) and hasattr(parsed_msg, dt_attr):
                                data_val = getattr(parsed_msg, df_attr)
                                data_type = getattr(parsed_msg, dt_attr)
                                
                                # pyubx2 should ideally handle signedness for I* types.
                                # For X3 (3-byte) fields like dataField in ESF-MEAS, it might be raw.
                                # If data_val is a positive integer representing a 24-bit two\'s complement negative number:
                                if data_type in [10,11,12,13,14,15] and (data_val & (1 << 23)): # Check 24th bit (sign bit for 24-bit int)
                                    data_val = data_val - (1 << 24) # Convert to signed from 24-bit representation

                                if data_type == 10:  # accX
                                    measurements["acc_x"] = data_val * (2**-10)  # Scale to m/s^2
                                    has_imu_data = True
                                elif data_type == 11:  # accY
                                    measurements["acc_y"] = data_val * (2**-10)
                                    has_imu_data = True
                                elif data_type == 12:  # accZ
                                    measurements["acc_z"] = data_val * (2**-10)
                                    has_imu_data = True
                                elif data_type == 13:  # gyroX
                                    measurements["gyro_x"] = data_val * (2**-12)  # Scale to rad/s
                                    has_imu_data = True
                                elif data_type == 14:  # gyroY
                                    measurements["gyro_y"] = data_val * (2**-12)
                                    has_imu_data = True
                                elif data_type == 15:  # gyroZ
                                    measurements["gyro_z"] = data_val * (2**-12)
                                    has_imu_data = True
                        
                        if has_imu_data:
                            measurements["time"] = current_time_for_msg
                            measurements["timeTag"] = parsed_msg.timeTag
                            esf_data_list.append(measurements)
                
                elif msg_type == "ESF-CAL":
                    current_time_for_msg = last_known_pvt_time
                    if current_time_for_msg:
                        sTtag = getattr(parsed_msg, 'sTtag', None)
                        cal_entry = {
                            "time": current_time_for_msg, # Approximate time from last NAV-PVT
                            "sTtag": sTtag
                        }

                        # REMOVED: esf_cal_type_map definition

                        items_found_in_message = 0
                        idx = 1 # Start with index _01 for dataField_XX and dataType_XX
                        while True:
                            df_attr = f"dataField_{idx:02d}"
                            dt_attr = f"dataType_{idx:02d}"

                            if hasattr(parsed_msg, dt_attr) and hasattr(parsed_msg, df_attr):
                                items_found_in_message +=1
                                data_type = getattr(parsed_msg, dt_attr)
                                raw_data_field = getattr(parsed_msg, df_attr) # Value from pyubx2 based on ESF-CAL definition
                                
                                field_name = f"cal_unknown_type_{data_type}" # Default if not found in map
                                value = raw_data_field # Default value is the raw_data_field

                                sensor_def = ESF_RAW_SENSOR_DEFINITION_MAP.get(data_type)

                                if sensor_def:
                                    field_name = sensor_def["name"]
                                    # Applying scaling from ESF_RAW_SENSOR_DEFINITION_MAP.
                                    # Note: raw_data_field for ESF-CAL is already processed by pyubx2,
                                    # not the 24-bit raw integer ESF_RAW_SENSOR_DEFINITION_MAP scales are designed for.
                                    # This direct application follows the user's request.
                                    value = raw_data_field * sensor_def["scale"]
                                # else: field_name and value remain as defaults for unknown types

                                cal_entry[field_name] = value
                                idx += 1
                            else:
                                # No more dataField_XX / dataType_XX attributes for this index
                                break 
                        
                        if items_found_in_message > 0:
                            esf_cal_list.append(cal_entry)

                elif msg_type == "ESF-INS":
                    current_time_for_msg = last_known_pvt_time
                    if current_time_for_msg and all(hasattr(parsed_msg, attr) for attr in [
                        'iTOW', 'xAngRate', 'yAngRate', 'zAngRate', 'xAccel', 'yAccel', 'zAccel',
                        'xAngRateValid', 'yAngRateValid', 'zAngRateValid', 
                        'xAccelValid', 'yAccelValid', 'zAccelValid'
                    ]):
                        # pyubx2 provides xAngRate in deg/s and xAccel in m/s^2
                        esf_ins_data = {
                            "time": current_time_for_msg,
                            "iTOW": parsed_msg.iTOW,
                            "ins_accel_x": parsed_msg.xAccel if parsed_msg.xAccelValid else np.nan,
                            "ins_accel_y": parsed_msg.yAccel if parsed_msg.yAccelValid else np.nan,
                            "ins_accel_z": parsed_msg.zAccel if parsed_msg.zAccelValid else np.nan,
                            "ins_gyro_x": parsed_msg.xAngRate if parsed_msg.xAngRateValid else np.nan,
                            "ins_gyro_y": parsed_msg.yAngRate if parsed_msg.yAngRateValid else np.nan,
                            "ins_gyro_z": parsed_msg.zAngRate if parsed_msg.zAngRateValid else np.nan
                        }
                        esf_ins_list.append(esf_ins_data)
                    elif not current_time_for_msg:
                        print(f"[DEBUG] Skipping ESF-INS (iTOW: {getattr(parsed_msg, 'iTOW', 'N/A')}) because last_known_pvt_time is None.")
                    # else: # Optionally handle cases where attributes are missing but time is known
                        # print(f"[DEBUG] Skipping ESF-INS (iTOW: {getattr(parsed_msg, 'iTOW', 'N/A')}) due to missing attributes.")

                elif msg_type == "ESF-ALG":
                    # print(f"Processing ESF-ALG message: {parsed_msg}") # Ensure this is active for debugging
                    current_time_for_msg = last_known_pvt_time
                    
                    required_attrs = ['yaw', 'pitch', 'roll', 'iTOW', 'status', 'tiltAlgError', 'yawAlgError']
                    all_attrs_present = all(hasattr(parsed_msg, attr) for attr in required_attrs)

                    if current_time_for_msg and all_attrs_present:
                        alg_status = parsed_msg.status
                        tilt_alg_error = parsed_msg.tiltAlgError
                        yaw_alg_error = parsed_msg.yawAlgError

                        esf_alg_list.append({
                            "time": current_time_for_msg,
                            "iTOW": parsed_msg.iTOW,
                            "yaw": parsed_msg.yaw,    # Degrees
                            "pitch": parsed_msg.pitch, # Degrees
                            "roll": parsed_msg.roll,   # Degrees
                            "algStatus": alg_status,
                            "tiltAlgError": tilt_alg_error,
                            "yawAlgError": yaw_alg_error,
                            # "angleError": parsed_msg.angleError, # Optionally store if needed for plotting/analysis
                            # "autoMntAlgOn": parsed_msg.autoMntAlgOn # Optionally store
                        })
                    else:
                        # Diagnostics for skipping
                        if not current_time_for_msg:
                            print(f"[DEBUG] Skipping ESF-ALG (iTOW: {getattr(parsed_msg, 'iTOW', 'N/A')}) because last_known_pvt_time is None.")
                        if not all_attrs_present:
                            missing_attrs = [attr for attr in required_attrs if not hasattr(parsed_msg, attr)]
                            print(f"[DEBUG] Skipping ESF-ALG (iTOW: {getattr(parsed_msg, 'iTOW', 'N/A')}) due to missing attributes: {missing_attrs}.")
                
                elif msg_type == "ESF-RAW":
                    current_time_for_msg = last_known_pvt_time # Approximate time of PVT
                    if current_time_for_msg:
                        data_entries = []
                        for attr_name in dir(parsed_msg):
                            if attr_name.startswith('data_') and attr_name[5:].isdigit():
                                data_entries.append(attr_name)
                        
                        grouped_raw_measurements = {}
                        
                        for data_attr in sorted(data_entries):
                            idx_str = data_attr.split('_')[1]
                            stag_attr = f"sTag_{idx_str}"
                            
                            if hasattr(parsed_msg, data_attr) and hasattr(parsed_msg, stag_attr):
                                data_bytes = getattr(parsed_msg, data_attr)
                                sensor_time_tag_bytes = getattr(parsed_msg, stag_attr)
                                
                                sensor_time_tag = int.from_bytes(sensor_time_tag_bytes, byteorder='little', signed=False)
                                
                                if len(data_bytes) >= 4:
                                    data_val_bytes = data_bytes[0:3]
                                    data_type_byte = data_bytes[3:4]
                                    
                                    raw_data_val_int = int.from_bytes(data_val_bytes, byteorder='little', signed=False)
                                    data_type = int.from_bytes(data_type_byte, byteorder='little', signed=False)
                                    
                                    # Handle 24-bit signed integer conversion
                                    if raw_data_val_int & (1 << 23):  # Check if sign bit (bit 23) is set
                                        signed_data_val = raw_data_val_int - (1 << 24)
                                    else:
                                        signed_data_val = raw_data_val_int
                                    
                                    sensor_def = ESF_RAW_SENSOR_DEFINITION_MAP.get(data_type)
                                    
                                    if sensor_def:
                                        field_name = sensor_def["name"]
                                        scale_factor = sensor_def["scale"]
                                        converted_val = signed_data_val * scale_factor
                                    else:
                                        field_name = f"unknown_type_0x{data_type:02x}"
                                        converted_val = signed_data_val # Store the signed raw value without scaling
                                    
                                    if sensor_time_tag not in grouped_raw_measurements:
                                        grouped_raw_measurements[sensor_time_tag] = {
                                            "time": current_time_for_msg,
                                            "sTag": sensor_time_tag,
                                        }
                                    
                                    grouped_raw_measurements[sensor_time_tag][field_name] = converted_val
                        
                        for measurement in grouped_raw_measurements.values():
                            esf_raw_stag_list.append(measurement)

        print(f"Finished parsing. Found {len(positions)} position, {len(satellite_data_list)} satellite, {len(esf_data_list)} ESF-MEAS, {len(esf_cal_list)} ESF-CAL, {len(esf_alg_list)} ESF-ALG, {len(esf_raw_list)} ESF-RAW records.")
        
        if esf_raw_stag_list:
            print("\\nSample ESF-RAW data (pre-interpolation, combined by sTag):")
            sample_display_list = esf_raw_stag_list[:min(2, len(esf_raw_stag_list))]

            for idx, sample_entry in enumerate(sample_display_list):
                print(f"--- Sample Entry {idx+1} ---")
                print(f"  Time (approx from PVT): {sample_entry['time']}")
                print(f"  sTag: {sample_entry['sTag']}")
                for key, value in sample_entry.items():
                    if key not in ['time', 'sTag']: 
                        print(f"  {key}: {value}")
        
        # Now that we've collected all ESF-RAW messages, perform timestamp interpolation
        # based on sTag values across the entire dataset
        if esf_raw_stag_list and first_pvt_time and last_pvt_time:
            print("\nInterpolating ESF-RAW timestamps based on sTag values...")
            # Find min and max sTag values
            min_stag = min([entry['sTag'] for entry in esf_raw_stag_list])
            max_stag = max([entry['sTag'] for entry in esf_raw_stag_list])
            stag_range = max_stag - min_stag
            
            if stag_range > 0:
                print(f"sTag range: {min_stag} to {max_stag} (range: {stag_range})")
                # Calculate the total time span of the dataset
                time_span = (last_pvt_time - first_pvt_time).total_seconds()
                print(f"Time span: {time_span:.3f} seconds")
                
                # Process each ESF-RAW entry and assign interpolated timestamp
                for entry in esf_raw_stag_list:
                    # Calculate the relative position of this sTag in the overall range
                    relative_pos = (entry['sTag'] - min_stag) / stag_range
                    
                    # Interpolate timestamp
                    time_offset = time_span * relative_pos
                    interpolated_time = first_pvt_time + timedelta(seconds=time_offset)
                    
                    # Update the timestamp
                    entry['time'] = interpolated_time
                    
                    # Add to the final list
                    esf_raw_list.append(entry)
                
                # Sort by timestamp to ensure proper sequence
                esf_raw_list.sort(key=lambda x: x['time'])
                
                print(f"Interpolated {len(esf_raw_list)} ESF-RAW timestamps")
                
                # Verify monotonically increasing timestamps
                timestamp_diffs = []
                for i in range(1, len(esf_raw_list)):
                    diff_ms = (esf_raw_list[i]['time'] - esf_raw_list[i-1]['time']).total_seconds() * 1000
                    timestamp_diffs.append(diff_ms)
                
                if all(diff >= 0 for diff in timestamp_diffs):
                    print("Verification: All timestamps are monotonically increasing ✓")
                else:
                    negative_diffs = sum(1 for diff in timestamp_diffs if diff < 0)
                    print(f"Warning: Found {negative_diffs} cases where timestamps go backwards")
                    
                # Calculate average interval
                if timestamp_diffs:
                    avg_interval = sum(timestamp_diffs) / len(timestamp_diffs)
                    print(f"Average interval between ESF-RAW timestamps: {avg_interval:.3f}ms")
            else:
                print("Error: Could not interpolate timestamps - all sTag values are identical or only one sTag value present")
                esf_raw_list = esf_raw_stag_list  # Use the collected data without interpolation
        elif esf_raw_stag_list: # sTag list exists, but PVT times might be missing or stag_range is 0
            if not first_pvt_time or not last_pvt_time:
                 print("Warning: ESF-RAW sTag data found, but PVT time references (first/last) are missing. Using unadjusted timestamps for ESF-RAW.")
            # If stag_range was 0, the message is already printed above.
            esf_raw_list = esf_raw_stag_list # Use as is, timestamps will be based on last_known_pvt_time at message arrival
        else:
            if not esf_raw_stag_list:
                print("No ESF-RAW messages with sTag to process for interpolation")

    except FileNotFoundError:
        print(f"Error: UBX file '{filename}' not found.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), raw_msgs
    except Exception as e:
        print(f"An error occurred during parsing: {e}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), raw_msgs
        
    df_positions = pd.DataFrame(positions)
    df_satellite = pd.DataFrame(satellite_data_list)
    df_esf = pd.DataFrame(esf_data_list)
    df_esf_cal = pd.DataFrame(esf_cal_list)
    df_esf_alg = pd.DataFrame(esf_alg_list)
    df_esf_raw = pd.DataFrame(esf_raw_list)
    df_esf_ins = pd.DataFrame(esf_ins_list) # Create DataFrame for ESF-INS
    df_nav_att = pd.DataFrame(nav_att_list) # Create DataFrame for NAV-ATT
    df_nav_dop = pd.DataFrame(nav_dop_list) # Create DataFrame for NAV-DOP

    return df_positions, df_satellite, df_esf, df_esf_cal, df_esf_alg, df_esf_raw, df_esf_ins, df_nav_att, df_nav_dop, raw_msgs

def _create_gps_1hz_csv(df_positions, df_nav_dop, output_dir):
    # Create GPS 1Hz DataFrame
    print("Creating GPS 1Hz DataFrame...")

    required_gps1_hz_cols = ['time', 'gnssFixOk', 'heading', 'speed', 'velN', 'velE', 'velD', 'lat', 'lon', 'hMSL', 'hAcc', 'vAcc', 'fixType']

    if not df_positions.empty:
        all(col in df_positions.columns for col in required_gps1_hz_cols)

        # Ensure time column is datetime object type
        df_positions['time'] = pd.to_datetime(df_positions['time'])

        # Downsample the data from 10 Hz to 1 Hz by taking every 10th row
        df_positions_1hz = df_positions[::10].copy()

        # Create the gps_1hz DataFrame
        gps_1hz = pd.DataFrame()
        
        # Convert timestamps to microseconds since the first timestamp
        if not df_positions_1hz['time'].empty:
            first_timestamp = df_positions_1hz['time'].iloc[0]
            gps_1hz['Timestamp'] = (df_positions_1hz['time'] - first_timestamp).dt.total_seconds() * 1_000_000
            gps_1hz['Timestamp'] = gps_1hz['Timestamp'].astype(int)
        else:
            gps_1hz['Timestamp'] = df_positions_1hz['time']
            print("Warning: Empty time column. No conversion applied.")
        
        # Add position data
        gps_1hz['Latitude'] = df_positions_1hz['lat']
        gps_1hz['Longitude'] = df_positions_1hz['lon']
        gps_1hz['Altitude'] = df_positions_1hz['hMSL']
        gps_1hz['Speed'] = df_positions_1hz['speed'] / 3.6  # km/h to m/s
        gps_1hz['Heading'] = df_positions_1hz['heading']
        
        # Add DOP data
        if not df_nav_dop.empty:
            # Ensure NAV-DOP DataFrame has the same time index
            df_nav_dop['time'] = pd.to_datetime(df_nav_dop['time'])
            gps_1hz = gps_1hz.merge(df_nav_dop, left_on='Timestamp', right_on='time', how='left')
            gps_1hz.drop(columns=['time'], inplace=True)
        else:
            gps_1hz['HDOP'] = 92.830002
            gps_1hz['VDOP'] = 26.480000
            gps_1hz['PDOP'] = 9651.000000
            gps_1hz['TDOP'] = 47.799999

        # Add accuracy estimates
        gps_1hz['hAcc'] = df_positions_1hz['hAcc']
        gps_1hz['vAcc'] = df_positions_1hz['vAcc']
        gps_1hz['sAcc'] = 0.311000

        # Add fix information
        gps_1hz['GPS fix'] = df_positions_1hz['fixType']
        gps_1hz['TTFF'] = 0.0
        gps_1hz['Fix valid'] = df_positions_1hz['gnssFixOk']
        
        # Save the navigation solution DataFrame to CSV
        gps_1hz.to_csv(os.path.join(output_dir, "gps_1hz.csv"), index=False, header=False)
        print("gps_1hz.csv created and saved.")
    else:
        print("Unable to create gps_1hz.csv: NAV-PVT data not available.")

def _create_imu104hz_csv(df_esf_raw, output_dir):
    # Create IMU 104Hz DataFrame
    print("Creating IMU 104Hz DataFrame...")
    required_imu_104hz_cols = ['time', 'raw_accel_x', 'raw_accel_y', 'raw_accel_z', 'raw_gyro_x', 'raw_gyro_y', 'raw_gyro_z']

    if not df_esf_raw.empty:
        all(col in df_esf_raw.columns for col in required_imu_104hz_cols)

        # Ensure time column is datetime object type
        df_esf_raw['time'] = pd.to_datetime(df_esf_raw['time'])

        # Create the imu_104hz DataFrame
        imu_104hz = pd.DataFrame()

        # Convert timestamps to microseconds since the first timestamp
        if not df_esf_raw['time'].empty:
            first_timestamp = df_esf_raw['time'].iloc[0]
            imu_104hz['Timestamp'] = (df_esf_raw['time'] - first_timestamp).dt.total_seconds() * 1_000_000
            imu_104hz['Timestamp'] = imu_104hz['Timestamp'].astype(int)
        else:
            imu_104hz['Timestamp'] = df_esf_raw['time']
            print("Warning: Empty time column. No conversion applied.")

        # Add IMU data
        imu_104hz['AccX'] = (df_esf_raw['raw_accel_x'] * 1000) / 9.81   # m/s^2 to milliG
        imu_104hz['AccY'] = (df_esf_raw['raw_accel_y'] * 1000) / 9.81   # m/s^2 to milliG
        imu_104hz['AccZ'] = (df_esf_raw['raw_accel_z'] * 1000) / 9.81   # m/s^2 to milliG
        imu_104hz['GyroX'] = df_esf_raw['raw_gyro_x'] * 1000    # DPS to milliDPS
        imu_104hz['GyroY'] = df_esf_raw['raw_gyro_y'] * 1000    # DPS to milliDPS
        imu_104hz['GyroZ'] = df_esf_raw['raw_gyro_z'] * 1000    # DPS to milliDPS

        # Convert columns to float
        imu_104hz['AccX'] = imu_104hz['AccX'].astype(np.float32)
        imu_104hz['AccY'] = imu_104hz['AccY'].astype(np.float32)
        imu_104hz['AccZ'] = imu_104hz['AccZ'].astype(np.float32)
        imu_104hz['GyroX'] = imu_104hz['GyroX'].astype(np.float32)
        imu_104hz['GyroY'] = imu_104hz['GyroY'].astype(np.float32)
        imu_104hz['GyroZ'] = imu_104hz['GyroZ'].astype(np.float32)
        

        # Save the IMU 104Hz DataFrame to CSV
        imu_104hz.to_csv(os.path.join(output_dir, "imu_104hz.csv"), index=False, header=False)
        print("imu_104hz.csv created and saved.")
    else:
        print("Unable to create imu_104hz.csv: ESF-RAW data not available.")

def convert_ubx_file_to_csv(ubx_file, output_dir):
    cache_file = ubx_file + ".pkl"

    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            df_positions, df_satellite, df_esf, df_esf_cal, df_esf_alg, df_esf_raw, df_esf_ins, df_nav_att, df_nav_dop, all_messages = pickle.load(f)
    else:
        df_positions, df_satellite, df_esf, df_esf_cal, df_esf_alg, df_esf_raw, df_esf_ins, df_nav_att, df_nav_dop, all_messages = parse_ubx_file(ubx_file)
        with open(cache_file, 'wb') as f:
            pickle.dump((df_positions, df_satellite, df_esf, df_esf_cal, df_esf_alg, df_esf_raw, df_esf_ins, df_nav_att, df_nav_dop, all_messages), f)

    # Convert to csv in test_csv/ folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    _create_gps_1hz_csv(df_positions, df_nav_dop, output_dir)
    _create_imu104hz_csv(df_esf_raw, output_dir)

    # Save all DataFrames to CSV files
    df_positions.to_csv(os.path.join(output_dir, "positions.csv"), index=False)
    df_satellite.to_csv(os.path.join(output_dir, "satellite.csv"), index=False)
    df_esf.to_csv(os.path.join(output_dir, "esf.csv"), index=False)
    df_esf_cal.to_csv(os.path.join(output_dir, "esf_cal.csv"), index=False)
    df_esf_alg.to_csv(os.path.join(output_dir, "esf_alg.csv"), index=False)
    df_esf_raw.to_csv(os.path.join(output_dir, "esf_raw.csv"), index=False)
    df_esf_ins.to_csv(os.path.join(output_dir, "esf_ins.csv"), index=False)
    df_nav_att.to_csv(os.path.join(output_dir, "nav_att.csv"), index=False)
    df_nav_dop.to_csv(os.path.join(output_dir, "nav_dop.csv"), index=False)

    print(f"All CSV files have been saved to the folder: {output_dir}/")

def plot_position_on_mapbox(df_pos, df_nav_att, mapbox_token): # Add df_nav_att
    if df_pos.empty:
        print("No position data to plot for Mapbox.")
        return
    if not all(col in df_pos.columns for col in ['lat', 'lon', 'time']):
        print("Position DataFrame missing 'lat', 'lon', or 'time' columns.")
        return
    
    if df_nav_att.empty or 'pitch' not in df_nav_att.columns or 'time' not in df_nav_att.columns:
        print("NAV-ATT DataFrame is empty or missing 'pitch' or 'time' columns. Plotting without pitch data.")
        # Fallback to fixType or a default color if NAV-ATT pitch is not available
        # For now, let's just use fixType as before if pitch is unavailable.
        fix_type_map = {
            0: "No Fix", 1: "Dead Reckoning", 2: "2D-Fix",
            3: "3D-Fix", 4: "GNSS + DR", 5: "Time Only"
        }
        df_pos['fixType_str'] = df_pos['fixType'].map(fix_type_map).fillna("Unknown")
        color_map = {
            "No Fix": "grey", "Dead Reckoning": "magenta", "2D-Fix": "yellow",
            "3D-Fix": "green", "GNSS + DR": "orange", "Time Only": "purple", "Unknown": "red"
        }
        fig = px.scatter_mapbox(df_pos,
                                lat="lat",
                                lon="lon",
                                color="fixType_str",
                                color_discrete_map=color_map,
                                hover_data=["time", "hMSL", "hAcc", "vAcc", "fixType_str"],
                                zoom=10,
                                height=600)
    else:
        print("Plotting position on Mapbox, colored by NAV-ATT pitch...")
        # Ensure 'time' columns are datetime objects for merging
        df_pos['time'] = pd.to_datetime(df_pos['time'])
        df_nav_att['time'] = pd.to_datetime(df_nav_att['time'])

        # Sort by time before asof merge
        df_pos = df_pos.sort_values('time')
        df_nav_att = df_nav_att.sort_values('time')

        # Merge df_pos with df_nav_att to get pitch at the closest time
        # Use asof merge to find the nearest NAV-ATT pitch for each position
        df_merged = pd.merge_asof(df_pos, 
                                  df_nav_att[['time', 'pitch']], 
                                  on='time', 
                                  direction='nearest', # or 'backward' or 'forward' depending on desired behavior
                                  tolerance=pd.Timedelta('1s')) # Optional: set a tolerance for matching

        if 'pitch' not in df_merged.columns or df_merged['pitch'].isnull().all():
            print("Could not merge pitch data effectively. Reverting to fixType coloring.")
            # Fallback to fixType coloring
            fix_type_map = {
                0: "No Fix", 1: "Dead Reckoning", 2: "2D-Fix",
                3: "3D-Fix", 4: "GNSS + DR", 5: "Time Only"
            }
            df_pos['fixType_str'] = df_pos['fixType'].map(fix_type_map).fillna("Unknown")
            color_map = {
                "No Fix": "grey", "Dead Reckoning": "magenta", "2D-Fix": "yellow",
                "3D-Fix": "green", "GNSS + DR": "orange", "Time Only": "purple", "Unknown": "red"
            }
            fig = px.scatter_mapbox(df_pos, # Use original df_pos if merge failed
                                    lat="lat",
                                    lon="lon",
                                    color="fixType_str",
                                    color_discrete_map=color_map,
                                    hover_data=["time", "hMSL", "hAcc", "vAcc", "fixType_str"],
                                    zoom=10,
                                    height=600)
        else:
            # Use pitch for color. Plotly will use a continuous colorscale.
            fig = px.scatter_mapbox(df_merged,
                                    lat="lat",
                                    lon="lon",
                                    color="pitch",  # Color by pitch
                                    color_continuous_scale=px.colors.sequential.Turbo, # Use Turbo colorscale
                                    range_color=[-5, 5], # Set range to -5 to 5 degrees
                                    hover_data=["time", "hMSL", "hAcc", "vAcc", "fixType", "pitch"], # Add pitch to hover
                                    zoom=10,
                                    height=600,
                                    title="Position Colored by NAV-ATT Pitch Angle (deg)")
            fig.update_layout(coloraxis_colorbar_title_text='Pitch (deg)')


    fig.update_layout(mapbox_accesstoken=mapbox_token, margin={"r":0,"t":0,"l":0,"b":0}, template="plotly_dark")
    config = {'scrollZoom': True}
    fig.show(config=config)

def plot_satellite_info(df_sat):
    if df_sat.empty:
        print("No satellite data to plot.")
        return

    print("Plotting satellite information...")

    cno_data_list = []
    if 'sats_details' in df_sat.columns and 'time' in df_sat.columns and 'iTOW' in df_sat.columns:
        for index, row in df_sat.iterrows():
            if isinstance(row['sats_details'], list):
                for sat_detail in row['sats_details']:
                    if isinstance(sat_detail, dict) and 'svId' in sat_detail and 'cno' in sat_detail:
                         cno_data_list.append({
                            'time': row['time'],
                            'iTOW': row['iTOW'],
                            'svId': sat_detail['svId'],
                            'cno': sat_detail['cno'],
                            'elev': sat_detail.get('elev', np.nan),
                            'azim': sat_detail.get('azim', np.nan)
                        })
    
    if cno_data_list:
        df_cno = pd.DataFrame(cno_data_list)
        if not df_cno.empty:
            df_cno['svId'] = df_cno['svId'].astype(str) 
            fig2 = px.scatter(df_cno, x="time", y="cno", color="svId", 
                              title="Satellite CNO vs Time (Approx.)",
                              hover_data=['svId', 'cno', 'elev', 'azim', 'iTOW'])
            fig2.update_layout(template="plotly_dark", xaxis_title="Approximate Time (from NAV-PVT)", yaxis_title="CNO (dBHz)")
            fig2.show()
    else:
        print("No detailed satellite CNO data to plot or relevant columns missing.")

raw_sensor_groups = [
    {
        "title": "Accelerometer (m/s²)",
        "components": ["raw_accel_x", "raw_accel_y", "raw_accel_z"],
        "names": ["Accel X", "Accel Y", "Accel Z"],
        "unit": "m/s²"
    },
    {
        "title": "Gyroscope (deg/s)",
        "components": ["raw_gyro_x", "raw_gyro_y", "raw_gyro_z"],
        "names": ["Gyro X", "Gyro Y", "Gyro Z"],
        "unit": "deg/s"
    },
    {
        "title": "Magnetometer (Raw)",
        "components": ["raw_mag_x", "raw_mag_y", "raw_mag_z"],
        "names": ["X", "Y", "Z"],
        "unit": "LSB" # Changed from uT to LSB to match ESF_RAW_SENSOR_DEFINITION_MAP
    },
    {
        "title": "Temperature (Raw)",
        "components": ["raw_temp"],
        "names": ["Temperature"],
        "unit": "°C"
    }
]

def plot_esf(df_esf_raw, df_esf_cal=None, df_positions=None, df_esf_alg=None, df_esf_ins=None, df_nav_att=None, df_ekf_attitude=None): # Add df_ekf_attitude parameter
    # Ensure raw_sensor_groups is accessible (e.g., global)
    # global raw_sensor_groups 

    subplot_configurations = [] # Stores configurations for each row

    # 1. ESF-RAW Sensor Plots (Accel, Gyro, Mag - excluding Temperature)
    if df_esf_raw is not None and not df_esf_raw.empty:
        for group_config in raw_sensor_groups:
            if "temperature" in group_config["title"].lower():
                continue 
            
            if all(col in df_esf_raw.columns and not df_esf_raw[col].dropna().empty for col in group_config["components"]):
                subplot_configurations.append({
                    "type": "standard_sensor",
                    "title": f"ESF-RAW {group_config['title']}", # Clarify source
                    "config": group_config,
                    "df": df_esf_raw
                })

    # 2. ESF-CAL Accelerometer Plot
    if df_esf_cal is not None and not df_esf_cal.empty and \
       all(col in df_esf_cal.columns and not df_esf_cal[col].dropna().empty for col in ['raw_accel_x', 'raw_accel_y', 'raw_accel_z']):
        subplot_configurations.append({
            "type": "standard_sensor",
            "title": "ESF-CAL Accelerometer (m/s²)",
            "config": { # group_info for this plot
                "title": "CAL Accel (m/s²)", # Y-axis title for this specific subplot
                "components": ["raw_accel_x", "raw_accel_y", "raw_accel_z"],
                "names": ["Accel X", "Accel Y", "Accel Z"],
                "unit": "m/s²"
            },
            "df": df_esf_cal
        })

    # 3. ESF-CAL Gyroscope Plot
    if df_esf_cal is not None and not df_esf_cal.empty and \
       all(col in df_esf_cal.columns and not df_esf_cal[col].dropna().empty for col in ['raw_gyro_x', 'raw_gyro_y', 'raw_gyro_z']):
        subplot_configurations.append({
            "type": "standard_sensor",
            "title": "ESF-CAL Gyroscope (deg/s)",
            "config": { # group_info for this plot
                "title": "CAL Gyro (deg/s)", # Y-axis title for this specific subplot
                "components": ["raw_gyro_x", "raw_gyro_y", "raw_gyro_z"],
                "names": ["Gyro X", "Gyro Y", "Gyro Z"],
                "unit": "deg/s"
            },
            "df": df_esf_cal
        })

    # 4. ESF-INS Accelerometer Plot
    if df_esf_ins is not None and not df_esf_ins.empty and \
       all(col in df_esf_ins.columns and not df_esf_ins[col].dropna().empty for col in ['ins_accel_x', 'ins_accel_y', 'ins_accel_z']):
        subplot_configurations.append({
            "type": "standard_sensor",
            "title": "ESF-INS Accelerometer (m/s²)", # Subplot title
            "config": { # group_info for this plot
                "title": "INS Accel (m/s²)", # Y-axis title for this specific subplot
                "components": ["ins_accel_x", "ins_accel_y", "ins_accel_z"],
                "names": ["Accel X", "Accel Y", "Accel Z"],
                "unit": "m/s²"
            },
            "df": df_esf_ins
        })

    # ESF-INS Gyroscope Plot
    if df_esf_ins is not None and not df_esf_ins.empty and \
       all(col in df_esf_ins.columns and not df_esf_ins[col].dropna().empty for col in ['ins_gyro_x', 'ins_gyro_y', 'ins_gyro_z']):
        subplot_configurations.append({
            "type": "standard_sensor",
            "title": "ESF-INS Gyroscope (deg/s)", # Subplot title
            "config": { # group_info for this plot
                "title": "INS Gyro (deg/s)", # Y-axis title for this specific subplot
                "components": ["ins_gyro_x", "ins_gyro_y", "ins_gyro_z"],
                "names": ["Gyro X", "Gyro Y", "Gyro Z"],
                "unit": "deg/s"
            },
            "df": df_esf_ins
        })

    # 2. Data for Combined Subplot ("Auxiliary Data")
    combined_plot_elements = []
    # Temperature from df_esf_raw
    temp_active = False
    if df_esf_raw is not None and 'raw_temp' in df_esf_raw.columns and not df_esf_raw['raw_temp'].dropna().empty:
        combined_plot_elements.append({
            "id": "temp", "name": "Temperature", "col": "raw_temp", "unit": "°C", 
            "df": df_esf_raw, "trace_type": "line", "yaxis_group": "temp"
        })
        temp_active = True

    # Horizontal Accuracy from df_positions
    hacc_active = False
    if df_positions is not None and 'hAcc' in df_positions.columns and not df_positions['hAcc'].dropna().empty:
        combined_plot_elements.append({
            "id": "hacc", "name": "hAcc", "col": "hAcc", "unit": "m", 
            "df": df_positions, "trace_type": "line", "yaxis_group": "hacc"
        })
        hacc_active = True

    # ALG Status from df_esf_alg
    alg_status_active = False
    if df_esf_alg is not None and 'algStatus' in df_esf_alg.columns and not df_esf_alg['algStatus'].dropna().empty:
        combined_plot_elements.append({
            "id": "alg_status", "name": "ALG Status", "col": "algStatus", "unit": "", 
            "df": df_esf_alg, "trace_type": "line", "yaxis_group": "fix_validity"
        })
        alg_status_active = True

    # ALG Errors from df_esf_alg
    alg_errors_active = False
    if df_esf_alg is not None and \
       ('tiltAlgError' in df_esf_alg.columns and not df_esf_alg['tiltAlgError'].dropna().empty) and \
       ('yawAlgError' in df_esf_alg.columns and not df_esf_alg['yawAlgError'].dropna().empty):
        combined_plot_elements.append({
            "id": "alg_tilt_error", "name": "Tilt Error", "col": "tiltAlgError", "unit": "deg", 
            "df": df_esf_alg, "trace_type": "line", "yaxis_group": "alg_error"
        })
        combined_plot_elements.append({
            "id": "alg_yaw_error", "name": "Yaw Error", "col": "yawAlgError", "unit": "deg", 
            "df": df_esf_alg, "trace_type": "line", "yaxis_group": "alg_error"
        })
        alg_errors_active = True

    # GNSS Fix OK from df_positions
    gnss_fix_ok_active = False
    if df_positions is not None and 'gnssFixOk' in df_positions.columns and not df_positions['gnssFixOk'].dropna().empty:
        combined_plot_elements.append({
            "id": "gnss_fix_ok", "name": "GNSS Fix OK", "col": "gnssFixOk", "unit": "", 
            "df": df_positions, "trace_type": "line", "yaxis_group": "fix_validity"
        })
        gnss_fix_ok_active = True
    
    if combined_plot_elements:
        subplot_configurations.append({
            "type": "combined_aux",
            "title": "Auxiliary Data", # Simplified title
            "elements": combined_plot_elements
        })

    # 3. ESF ALG Angles Plot
    alg_angles_active = False
    if df_esf_alg is not None and not df_esf_alg.empty and \
       all(col in df_esf_alg.columns and not df_esf_alg[col].dropna().empty for col in ['roll', 'pitch', 'yaw']):
        subplot_configurations.append({
            "type": "alg_angles",
            "title": "ESF ALG Angles (deg)",
            "df": df_esf_alg
        })
        alg_angles_active = True

    # /* START: Comment out old combined NAV-ATT & EKF plot configuration
    # # 4. NAV-ATT Angles Plot
    # nav_att_angles_active = False
    # if df_nav_att is not None and not df_nav_att.empty and \\\n       all(col in df_nav_att.columns and not df_nav_att[col].dropna().empty for col in [\'roll\', \'pitch\', \'heading\']):
    #     subplot_configurations.append({
    #         "type": "nav_att_angles",
    #         "title": "NAV-ATT & EKF Angles (deg)", # Updated title
    #         "df_nav_att": df_nav_att, # Pass NAV-ATT data
    #         "df_ekf_attitude": df_ekf_attitude # Pass EKF data
    #     })
    #     nav_att_angles_active = True # Keep this for logic if needed, though plot presence is now conditional inside
    # elif df_ekf_attitude is not None and not df_ekf_attitude.empty and \\\n         all(col in df_ekf_attitude.columns and not df_ekf_attitude[col].dropna().empty for col in [\'ekf_roll\', \'ekf_pitch\', \'ekf_yaw\']):
    #     # If only EKF data is available, still create the subplot for it
    #     subplot_configurations.append({
    #         "type": "nav_att_angles", # Reuse the type, logic will handle missing NAV-ATT
    #         "title": "EKF Angles (deg)", # Title for EKF only
    #         "df_nav_att": None, # NAV-ATT is not available
    #         "df_ekf_attitude": df_ekf_attitude # Pass EKF data
    #     })
    #     # nav_att_angles_active can remain false or be set true; plot logic handles data presence
    # END: Comment out old combined NAV-ATT & EKF plot configuration */

    # New individual angle plots
    nav_att_available = df_nav_att is not None and not df_nav_att.empty
    ekf_attitude_available = df_ekf_attitude is not None and not df_ekf_attitude.empty

    # Roll Plot
    # Check if there's any roll data to plot (either NAV or EKF)
    can_plot_roll = False
    if nav_att_available and 'roll' in df_nav_att.columns and not df_nav_att['roll'].dropna().empty:
        can_plot_roll = True
    elif ekf_attitude_available and 'ekf_roll' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_roll'].dropna().empty:
        can_plot_roll = True
    
    if can_plot_roll:
        subplot_configurations.append({
            "type": "roll_angle_plot",
            "title": "Roll Angle (deg)",
            "df_nav_att": df_nav_att if nav_att_available and 'roll' in df_nav_att.columns and not df_nav_att['roll'].dropna().empty else None,
            "df_ekf_attitude": df_ekf_attitude if ekf_attitude_available and 'ekf_roll' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_roll'].dropna().empty else None
        })

    # Pitch Plot
    can_plot_pitch = False
    if nav_att_available and 'pitch' in df_nav_att.columns and not df_nav_att['pitch'].dropna().empty:
        can_plot_pitch = True
    elif ekf_attitude_available and 'ekf_pitch' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_pitch'].dropna().empty:
        can_plot_pitch = True

    if can_plot_pitch:
        subplot_configurations.append({
            "type": "pitch_angle_plot",
            "title": "Pitch Angle (deg)",
            "df_nav_att": df_nav_att if nav_att_available and 'pitch' in df_nav_att.columns and not df_nav_att['pitch'].dropna().empty else None,
            "df_ekf_attitude": df_ekf_attitude if ekf_attitude_available and 'ekf_pitch' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_pitch'].dropna().empty else None
        })

    # Yaw/Heading Plot
    can_plot_yaw = False
    if nav_att_available and 'heading' in df_nav_att.columns and not df_nav_att['heading'].dropna().empty:
        can_plot_yaw = True
    elif ekf_attitude_available and 'ekf_yaw' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_yaw'].dropna().empty:
        can_plot_yaw = True
        
    if can_plot_yaw:
        subplot_configurations.append({
            "type": "yaw_angle_plot",
            "title": "Yaw/Heading Angle (deg)",
            "df_nav_att": df_nav_att if nav_att_available and 'heading' in df_nav_att.columns and not df_nav_att['heading'].dropna().empty else None,
            "df_ekf_attitude": df_ekf_attitude if ekf_attitude_available and 'ekf_yaw' in df_ekf_attitude.columns and not df_ekf_attitude['ekf_yaw'].dropna().empty else None
        })

    # 5. Ground Speed Plot from df_positions
    if df_positions is not None and not df_positions.empty and \
       'speed' in df_positions.columns and not df_positions['speed'].dropna().empty:
        subplot_configurations.append({
            "type": "standard_sensor", # Reusing standard_sensor type
            "title": "Ground Speed (km/h)", # Subplot title
            "config": {
                "title": "Speed (km/h)", # Y-axis title
                "components": ["speed"],
                "names": ["Ground Speed"],
                "unit": "km/h"
            },
            "df": df_positions
        })

    if not subplot_configurations:
        print("No data to plot in plot_esf.")
        return

    num_cols = 2
    total_rows = math.ceil(len(subplot_configurations) / num_cols)
    subplot_titles = [config["title"] for config in subplot_configurations]
    
    # Adjust row heights (example, can be fine-tuned)
    row_heights = [1.0 / total_rows] * total_rows if total_rows > 0 else None

    fig = make_subplots(
        rows=total_rows,
        cols=num_cols, # Use 2 columns
        shared_xaxes=True,
        subplot_titles=subplot_titles,
        row_heights=row_heights,
        vertical_spacing=0.05 if total_rows > 1 else 0.02,
        horizontal_spacing=0.05 # Add horizontal spacing
    )

    layout_yaxes_config = {} # To store all y-axis configurations for fig.update_layout
    next_global_secondary_yaxis_idx = (total_rows * num_cols) + 1 # Adjust for multiple columns

    for i, plot_config in enumerate(subplot_configurations):
        current_row_idx = math.ceil((i + 1) / num_cols) # Calculate row based on 2 columns
        current_col_idx = (i % num_cols) + 1 # Calculate col based on 2 columns
        
        yaxis_idx_for_cell = (current_col_idx - 1) * total_rows + current_row_idx
        primary_yaxis_name_for_row = f'yaxis{yaxis_idx_for_cell if yaxis_idx_for_cell > 1 else ""}'

        if plot_config["type"] == "standard_sensor":
            group_info = plot_config["config"]
            df_current = plot_config["df"]
            for comp_idx, component in enumerate(group_info["components"]):
                fig.add_trace(
                    go.Scatter(x=df_current['time'], y=df_current[component], mode='lines', 
                               name=f'{group_info["names"][comp_idx]} ({group_info["unit"]})'),
                    row=current_row_idx, col=current_col_idx
                )
            # Configure the primary y-axis for this standard sensor plot
            layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}

        elif plot_config["type"] == "combined_aux":
            elements = plot_config["elements"]
            # For the combined plot, assign traces to y, y2, y3... within its row
            # Map a yaxis_group (e.g., "temp", "hacc", "alg_error") to a cell-local y-axis specifier ('y', 'y2', 'y3')
            cell_local_yaxis_assignment = {}
            current_cell_local_yaxis_idx = 1 # 1 for 'y', 2 for 'y2', etc.

            for el_idx, element_config in enumerate(elements):
                el_yaxis_group = element_config["yaxis_group"]
                trace_yaxis_specifier_for_cell = '' # Default for primary y-axis of the cell

                if el_yaxis_group not in cell_local_yaxis_assignment:
                    if current_cell_local_yaxis_idx == 1:
                        trace_yaxis_specifier_for_cell = 'y' # Primary for this cell
                    else:
                        trace_yaxis_specifier_for_cell = f'y{current_cell_local_yaxis_idx}' # y2, y3...
                    cell_local_yaxis_assignment[el_yaxis_group] = trace_yaxis_specifier_for_cell
                    current_cell_local_yaxis_idx += 1
                else:
                    trace_yaxis_specifier_for_cell = cell_local_yaxis_assignment[el_yaxis_group]
                
                df_el = element_config["df"]
                el_name = element_config["name"]
                el_col = element_config["col"]
                el_unit = element_config["unit"]

                fig.add_trace(
                    go.Scatter(x=df_el['time'], y=df_el[el_col], mode='lines', 
                               name=f"{el_name} ({el_unit})" if el_unit else el_name, 
                               yaxis=trace_yaxis_specifier_for_cell),
                    row=current_row_idx, col=current_col_idx
                )

                # Configure the global y-axis this trace maps to
                y_title_text = f"{el_name.split(' ')[0]} ({el_unit})" if el_unit else el_name.split(' ')[0]
                if el_yaxis_group == "alg_error" and "Tilt Error" in el_name : # Title for the group
                     y_title_text = f"ALG Error ({el_unit})"
                elif el_yaxis_group == "alg_error" and "Yaw Error" in el_name:
                    continue # Title already set by Tilt Error for the group
                elif el_yaxis_group == "fix_validity":
                    y_title_text = "Fix Valid"

                if trace_yaxis_specifier_for_cell == 'y': # Maps to the row's primary y-axis
                    if primary_yaxis_name_for_row not in layout_yaxes_config: # Set title if not already set by another trace on same primary
                        layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}
                else: # Maps to a secondary y-axis for this row (y2, y3...)
                    # This needs a new global y-axis definition in the layout
                    global_secondary_yaxis_name = f"yaxis{next_global_secondary_yaxis_idx}"
                    if global_secondary_yaxis_name not in layout_yaxes_config:
                        layout_yaxes_config[global_secondary_yaxis_name] = {
                            "title": "",
                            "side": "right",
                            "overlaying": primary_yaxis_name_for_row.replace('yaxis','y'), # Overlay on primary of this cell
                            "anchor": f"x{yaxis_idx_for_cell if yaxis_idx_for_cell > 1 else ''}", # Anchor to the correct x-axis of the cell
                            "autoshift": True,
                            "showgrid": False # Often good for secondary axes
                        }
                        next_global_secondary_yaxis_idx += 1
            
        elif plot_config["type"] == "alg_angles":
            df_current = plot_config["df"]
            angle_components = [("roll", "Roll"), ("pitch", "Pitch"), ("yaw", "Yaw")]
            for comp_key, comp_name in angle_components:
                fig.add_trace(
                    go.Scatter(x=df_current['time'], y=df_current[comp_key], mode='lines', name=f'{comp_name} (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}
            
        # /* START: Comment out old combined NAV-ATT & EKF plot rendering
        # elif plot_config["type"] == "nav_att_angles":
        #     df_nav = plot_config.get("df_nav_att")
        #     df_ekf = plot_config.get("df_ekf_attitude")
        #     plot_title_set = False

        #     if df_nav is not None and not df_nav.empty:
        #         angle_components_nav = [("roll", "NAV Roll"), ("pitch", "NAV Pitch"), ("heading", "NAV Heading")]
        #         for comp_key, comp_name in angle_components_nav:
        #             if comp_key in df_nav.columns:
        #                 fig.add_trace(
        #                     go.Scatter(x=df_nav['time'], y=df_nav[comp_key], mode='lines', name=f'{comp_name} (deg)'),
        #                     row=current_row_idx, col=current_col_idx
        #                 )
        #         if not plot_title_set:
        #             layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "Angle (deg)", "side": "left"}
        #             plot_title_set = True
            
        #     if df_ekf is not None and not df_ekf.empty:
        #         angle_components_ekf = [("ekf_roll", "EKF Roll"), ("ekf_pitch", "EKF Pitch"), ("ekf_yaw", "EKF Yaw")]
        #         for idx, (comp_key, comp_name) in enumerate(angle_components_ekf):
        #             if comp_key in df_ekf.columns:
        #                 fig.add_trace(
        #                     go.Scatter(x=df_ekf['time'], y=df_ekf[comp_key], mode='lines', name=f'{comp_name} (deg)'),
        #                     row=current_row_idx, col=current_col_idx
        #                 )
        #         if not plot_title_set:
        #             layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "Angle (deg)", "side": "left"}
        #             plot_title_set = True
            
        #     if not plot_title_set:
        #          layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "Angle (deg)", "side": "left"}
        # END: Comment out old combined NAV-ATT & EKF plot rendering */

        elif plot_config["type"] == "roll_angle_plot":
            df_nav = plot_config.get("df_nav_att")
            df_ekf = plot_config.get("df_ekf_attitude")
            
            if df_nav is not None and 'roll' in df_nav.columns: # Already checked for non-empty in config
                fig.add_trace(
                    go.Scatter(x=df_nav['time'], y=df_nav['roll'], mode='lines', name='NAV Roll (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            if df_ekf is not None and 'ekf_roll' in df_ekf.columns: # Already checked for non-empty in config
                fig.add_trace(
                    go.Scatter(x=df_ekf['time'], y=df_ekf['ekf_roll'], mode='lines', name='EKF Roll (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}

        elif plot_config["type"] == "pitch_angle_plot":
            df_nav = plot_config.get("df_nav_att")
            df_ekf = plot_config.get("df_ekf_attitude")

            if df_nav is not None and 'pitch' in df_nav.columns:
                fig.add_trace(
                    go.Scatter(x=df_nav['time'], y=df_nav['pitch'], mode='lines', name='NAV Pitch (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            if df_ekf is not None and 'ekf_pitch' in df_ekf.columns:
                fig.add_trace(
                    go.Scatter(x=df_ekf['time'], y=df_ekf['ekf_pitch'], mode='lines', name='EKF Pitch (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}

        elif plot_config["type"] == "yaw_angle_plot":
            df_nav = plot_config.get("df_nav_att")
            df_ekf = plot_config.get("df_ekf_attitude")

            if df_nav is not None and 'heading' in df_nav.columns:
                fig.add_trace(
                    go.Scatter(x=df_nav['time'], y=df_nav['heading'], mode='lines', name='NAV Heading (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            if df_ekf is not None and 'ekf_yaw' in df_ekf.columns:
                fig.add_trace(
                    go.Scatter(x=df_ekf['time'], y=df_ekf['ekf_yaw'], mode='lines', name='EKF Yaw (deg)'),
                    row=current_row_idx, col=current_col_idx
                )
            layout_yaxes_config[primary_yaxis_name_for_row] = {"title": "", "side": "left"}


    fig.update_layout(**layout_yaxes_config)
    fig.update_layout(
        height=max(600, (280 * total_rows) + 50 * (total_rows -1) ),
        title_text="ESF Raw Sensor Data, Auxiliary Info & ALG Angles", 
        template="plotly_dark",
        legend_orientation="h", legend_yanchor="bottom", legend_y=1.07, legend_xanchor="right", legend_x=1,
        xaxis2_matches='x1',
        xaxis4_matches='x1',
        xaxis6_matches='x1',
        xaxis8_matches='x1',
        xaxis10_matches='x1',
        margin=dict(t=400) # Add top margin to prevent overlap
    )
    fig.show()